#include <string.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpRender.h>
#if USE_FFMPEG
  #include <avformat.h>
  #include <swscale.h>

// add a video output stream
static AVStream *add_video_stream(AVFormatContext *oc, enum CodecID codec_id, int width, int height, int rate);
static AVStream *add_video_stream(AVFormatContext *oc, enum CodecID codec_id, int width, int height, int rate){
  AVStream       *st;
  AVCodecContext *c;

  st=av_new_stream(oc,0);
  if(!st)
    SID_trap_error("Could not alloc stream",ERROR_LOGIC);
  c               =st->codec;
  c->codec_id     =codec_id;
  c->codec_type   =CODEC_TYPE_VIDEO;
  c->width        =width+width%2;
  c->height       =height+height%2;
  c->bit_rate     =16*1024*1024;
/*
  c->qscale       =MAX(1,MIN(31,GBPGFX_QSCALE));
  c->chroma_qscale=c->chroma_qscale_table[c->qscale];
  c->y_dc_scale   =c->y_dc_scale_table[c->qscale];
  c->c_dc_scale   =c->c_dc_scale_table[c->chroma_qscale]; 
*/

  // time base: this is the fundamental unit of time (in seconds) in terms
  //   of which frame timestamps are represented. for fixed-fps content,
  //   timebase should be 1/framerate and timestamp increments should be
  //   identically 1.
  c->time_base.den=rate;
  c->time_base.num=1;
  c->gop_size     =GBPGFX_GOP_SIZE; // Emit one intra frame every GBPGFX_GOP_SIZE frames at most
  c->pix_fmt      =PIX_FMT_YUV420P;

  if (c->codec_id==CODEC_ID_MPEG2VIDEO){
    // just for testing, we also add B frames
    c->max_b_frames=2;
  }

  if(c->codec_id==CODEC_ID_MPEG1VIDEO){
    // needed to avoid using macroblocks in which some coeffs overflow
    //   this doesnt happen with normal video, it just happens here as the
    //   motion of the chroma plane doesnt match the luma plane
    c->mb_decision=2;
  }

  // some formats want stream headers to be separate
  if(!strcmp(oc->oformat->name, "mp4") || !strcmp(oc->oformat->name, "mov") || !strcmp(oc->oformat->name, "3gp"))
    c->flags |= CODEC_FLAG_GLOBAL_HEADER;

  return st;
}

static AVFrame *alloc_picture(int pix_fmt, int width, int height);
static AVFrame *alloc_picture(int pix_fmt, int width, int height){
  AVFrame *picture;
  uint8_t *picture_buf;
  int      size;
  picture = avcodec_alloc_frame();
  if(!picture)
    return NULL;
  size       =avpicture_get_size(pix_fmt, width, height);
  picture_buf=(uint8_t *)av_malloc(size);
  if(!picture_buf) {
    av_free(picture);
    return NULL;
  }
  avpicture_fill((AVPicture *)picture, picture_buf,
                 pix_fmt, width, height);
  return picture;
}

static void open_video(movie_info *movie);
static void open_video(movie_info *movie){
  AVCodecContext  *c;
  AVFormatContext *oc;
  AVStream        *st;
  AVCodec         *codec;

  c =movie->video_stream->codec;
  oc=movie->video_context;
  st=movie->video_stream;

  // Find the video encoder
  codec=avcodec_find_encoder(c->codec_id);
  if(!codec)
    SID_trap_error("Codec not found",ERROR_LOGIC);

  // Open the codec
  if(avcodec_open(c,codec)<0)
    SID_trap_error("Could not open codec",ERROR_LOGIC);

  movie->video_outbuf     =NULL;
  movie->video_outbuf_size=0;
  if (!(oc->oformat->flags & AVFMT_RAWPICTURE)) {
    // allocate output buffer 
    // XXX: API change will be done 
    // buffers passed into lav* can be allocated any way you prefer,
    //   as long as they're aligned enough for the architecture, and
    //   they're freed appropriately (such as using av_free for buffers
    //   allocated with av_malloc)
    movie->video_outbuf_size=256*SIZE_OF_KILOBYTE;
    movie->video_outbuf     =(uint8_t *)av_malloc(movie->video_outbuf_size);
  }

  // Allocate the encoded raw picture
  movie->picture=alloc_picture(c->pix_fmt,c->width,c->height);
  if(!movie->picture)
    SID_trap_error("Could not allocate picture",ERROR_LOGIC);

  // if the output format is not YUV420P, then a temporary YUV420P
  //   picture is needed too. It is then converted to the required
  //   output format
  movie->temp_picture=NULL;
  if(c->pix_fmt!=PIX_FMT_YUV420P) {
    movie->temp_picture=alloc_picture(PIX_FMT_YUV420P,c->width,c->height);
    if (!movie->temp_picture) 
      SID_trap_error("Could not allocate temporary picture",ERROR_LOGIC);
  }
}
#endif

void open_movie(char       *filename,
                int         width,
                int         height,
                int         n_frames,
                int         frame_rate,
                movie_info *movie){
#if USE_FFMPEG
  AVOutputFormat *fmt;
  
  if(SID.I_am_Master){
    SID_log("Initializing movie {%s}...",SID_LOG_OPEN,filename);

    // Quiet all those pesky av_log messages
    av_log_set_level(GBPGFX_AV_LOG_LEVEL);
    SID_log("(av_log_level set to %d)...",SID_LOG_CONTINUE,av_log_get_level());

    // Set frame rate, etc.
    movie->frame_count   =0;
    movie->n_frames      =n_frames;
    movie->frame_rate    =frame_rate;
    movie->total_duration=(double)n_frames/(double)frame_rate;

    // Initialize libavcodec, and register all codecs and formats
    av_register_all();

    // Auto detect the output format from the name. Default is mpeg.
    fmt=guess_format(NULL, filename, NULL);
    if(!fmt){
      SID_log("file extension unknown; using default=%s...",SID_LOG_CONTINUE,GBPGFX_FORMAT_DEFAULT);
      fmt = guess_format(GBPGFX_FORMAT_DEFAULT, NULL, NULL);
      if(!fmt)
        SID_trap_error("Default format is invalid",ERROR_LOGIC);
    }
    if(!fmt)
      SID_trap_error("Could not find suitable output format",ERROR_LOGIC);

    // Allocate the output media context
    movie->video_context=avformat_alloc_context();
    if(!movie->video_context)
      SID_trap_error("Movie context allocation error",ERROR_MEMORY);
    movie->video_context->oformat=fmt;
    snprintf(movie->video_context->filename,sizeof(movie->video_context->filename),"%s",filename);

    // Initialize the codec
    movie->video_stream=NULL;
    if(fmt->video_codec!=CODEC_ID_NONE)
      movie->video_stream=add_video_stream(movie->video_context,fmt->video_codec,width,height,movie->frame_rate);

    // Set the output parameters (must be done even if no parameters)
    if(av_set_parameters(movie->video_context,NULL)<0)
      SID_trap_error("Invalid format parameters",ERROR_LOGIC);
    dump_format(movie->video_context,0,filename,1);

    // Open the video codec and allocate the necessary encode buffers
    if(movie->video_stream)
      open_video(movie);

    // Open the output file, if needed
    if(!(fmt->flags & AVFMT_NOFILE)) {
      if(url_fopen(&movie->video_context->pb,filename,URL_WRONLY)<0)
        SID_trap_error("Could not open file {%s}",ERROR_IO_OPEN,filename);
    }

    // Write the stream header, if any
    av_write_header(movie->video_context);

    SID_log("Done.",SID_LOG_CLOSE);
  }
#else
  SID_trap_error("Function not supported.  FFMPEG is not installed.",ERROR_LOGIC);
#endif
}
