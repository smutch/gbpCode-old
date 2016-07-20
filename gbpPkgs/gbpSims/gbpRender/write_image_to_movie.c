#include <gd.h>
#include <gbpLib.h>
#include <gbpRender.h>

#if USE_FFMPEG
#include <avformat.h>
#include <swscale.h>

// Prepare the frames
static void fill_yuv_image(AVFrame *picture, image_info *image, AVCodecContext *c);
static void fill_yuv_image(AVFrame *picture, image_info *image, AVCodecContext *c){
  int           i_x;
  int           i_y;
  int           i;
  int          *dest;
  AVFrame      *tmp_frame;
  static struct SwsContext *img_convert_ctx;

  tmp_frame = avcodec_alloc_frame();
  avpicture_alloc((AVPicture*)tmp_frame, PIX_FMT_RGBA32, c->width, c->height);

  // Convert gdlib image to the FFmpeg image buffer
  dest=(int*)tmp_frame->data[0];
  for(i_y=0;i_y<c->height;i_y++) {
    for(i_x=0;i_x<c->width;i_x++) {
      dest[i_x]=image->gd_ptr->tpixels[i_y][i_x];
    }
    dest+=c->width;
  }

  // convert to the needed format
  img_convert_ctx = sws_getContext(c->width, c->height,
                                   PIX_FMT_RGBA32,
                                   c->width, c->height,
                                   c->pix_fmt,
                                   GBPGFX_SWS_FLAGS, NULL, NULL, NULL);
  sws_scale(img_convert_ctx, 
            tmp_frame->data, 
            tmp_frame->linesize, 
            0, 
            c->height, 
            picture->data, 
            picture->linesize);

}
#endif

void write_image_to_movie(image_info *image, movie_info *movie){
#if USE_FFMPEG
  AVCodecContext *c;
  int             out_size;
  int             ret;
  static struct   SwsContext *img_convert_ctx;
  int             flag_silent;

  if(SID.I_am_Master){

    c = movie->video_stream->codec;

    flag_silent=FALSE;
    if(movie->frame_count>=movie->n_frames) {
      // no more frames to compress. The codec has a latency of a few
      //   frames if using B frames, so we get the last frames by
      //   passing the same picture again 
      flag_silent=TRUE;
    } 
    else {
      SID_log("Writing frame %d of %d to movie...",SID_LOG_OPEN,movie->frame_count+1,movie->n_frames);
      if(c->pix_fmt!=PIX_FMT_YUV420P) {
        // as we only generate a YUV420 picture, we must convert it to the codec pixel format if needed
        if (img_convert_ctx == NULL) {
          img_convert_ctx = sws_getContext(c->width, c->height,
                                           PIX_FMT_YUV420P,
                                           c->width, c->height,
                                           c->pix_fmt,
                                           GBPGFX_SWS_FLAGS, NULL, NULL, NULL);
          fprintf(stderr,"test1\n");
          if(img_convert_ctx==NULL)
            SID_trap_error("Cannot initialize the conversion context",ERROR_LOGIC);
        }
        fill_yuv_image(movie->temp_picture,image,c);
        sws_scale(img_convert_ctx,
                  movie->temp_picture->data,
                  movie->temp_picture->linesize,
                  0, 
                  c->height,
                  movie->picture->data,
                  movie->picture->linesize);
      } 
      else
        fill_yuv_image(movie->picture,image,c);
    }

    // Raw video case. The API will change slightly in the near future for that 
    if(movie->video_context->oformat->flags & AVFMT_RAWPICTURE) {
      AVPacket pkt;
      av_init_packet(&pkt);
      pkt.flags       |= PKT_FLAG_KEY;
      pkt.stream_index = movie->video_stream->index;
      pkt.data         = (uint8_t *)movie->picture;
      pkt.size         = sizeof(AVPicture);
      ret              = av_write_frame(movie->video_context,&pkt);
    } 
    // encode the image
    else {
      out_size=avcodec_encode_video(c,movie->video_outbuf,movie->video_outbuf_size,movie->picture);
      // if zero size, it means the image was buffered
      if(out_size>0) {
        AVPacket pkt;
        av_init_packet(&pkt);
        pkt.pts= av_rescale_q(c->coded_frame->pts, c->time_base, movie->video_stream->time_base);
        if(c->coded_frame->key_frame)
          pkt.flags |= PKT_FLAG_KEY;
        pkt.stream_index=movie->video_stream->index;
        pkt.data        =movie->video_outbuf;
        pkt.size        =out_size;

        // write the compressed frame in the media file
        ret = av_write_frame(movie->video_context,&pkt);
      } 
      else
        ret=0;
    }
    if(ret!=0)
      SID_trap_error("Error while writing video frame",ERROR_LOGIC);
    movie->frame_count++;
    if(!flag_silent)
      SID_log("Done.",SID_LOG_CLOSE);
  }
#else
  SID_trap_error("Routine not supported.  FFMPEG not installed.",ERROR_LOGIC);
#endif
}
