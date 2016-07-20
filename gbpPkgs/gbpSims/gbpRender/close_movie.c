#include <gbpLib.h>
#include <gbpRender.h>
#if USE_FFMPEG
  #include "libavformat/avformat.h"
  #include "libswscale/swscale.h"
#endif

void close_movie(movie_info *movie){
  int i;
  #if USE_FFMPEG
  if(SID.I_am_Master){
    // The codec can have a latency of a few
    //   frames if using B frames, so we get the last frames by
    //   passing the same picture again
    /*
    while((double)movie->video_stream->pts.val*movie->video_stream->time_base.num/movie->video_stream->time_base.den<movie->total_duration)
      write_image_to_movie(NULL,movie);
    */

    // Close the codec and free image buffers
    if(movie->video_stream){
      avcodec_close(movie->video_stream->codec);
      av_free(movie->picture->data[0]);
      av_free(movie->picture);
      if(movie->temp_picture){
        av_free(movie->temp_picture->data[0]);
        av_free(movie->temp_picture);
      }
      av_free(movie->video_outbuf);
    }

    // Write the trailer, if any
    av_write_trailer(movie->video_context);

    // Free all streams
    for(i = 0; i < movie->video_context->nb_streams; i++) {
      av_freep(&movie->video_context->streams[i]->codec);
      av_freep(&movie->video_context->streams[i]);
    }

    // Close the output file
    if(!(movie->video_context->oformat->flags & AVFMT_NOFILE))
      url_fclose(movie->video_context->pb);

    // Free the context
    av_free(movie->video_context);
  
    movie->video_stream     =NULL;
    movie->video_context    =NULL;
    movie->picture          =NULL;
    movie->temp_picture     =NULL;
    movie->video_outbuf     =NULL;
    movie->video_outbuf_size=0;
    movie->width            =0;
    movie->height           =0;
    movie->n_pixels         =0;
    movie->n_frames         =0;
    movie->frame_rate       =0;
  }
  #else
    SID_trap_error("Routine not supported.  FFMPEG not installed.",ERROR_LOGIC);
  #endif
}
