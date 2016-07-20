#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
int main(int argc, char *argv[]){

    // Parse cmd line
    if(argc<9)
       SID_trap_error("Invalid syntax.",ERROR_SYNTAX);
    SID_init(&argc,&argv,NULL,NULL);

    char   filename[MAX_FILENAME_LENGTH];
    sprintf(filename,  argv[1]);
    int    x_col =atoi(argv[2]);
    int    y_col =atoi(argv[3]);
    int    z_col =atoi(argv[4]);
    float  x_cen =atof(argv[5]);
    float  y_cen =atof(argv[6]);
    float  z_cen =atof(argv[7]);
    float  radius=atof(argv[8]);
    float  radius2=radius*radius;

    // Process file
    char   *line       =NULL;
    size_t  line_length=0;
    FILE   *fp_in      =fopen(filename,"r");
    int     n_lines    =count_lines_data(fp_in);
    for(int i_line=0;i_line<n_lines;i_line++){
       float x_i;
       float y_i;
       float z_i;
       grab_next_line_data(fp_in,&line,&line_length);
       grab_float(line,x_col,&x_i);
       grab_float(line,y_col,&y_i);
       grab_float(line,z_col,&z_i);
       x_i-=x_cen;
       y_i-=y_cen;
       z_i-=z_cen;
       float d2=x_i*x_i+y_i*y_i+z_i*z_i;
       if(d2<radius2) printf("%s",line);
    }
    fclose(fp_in);

    SID_exit(ERROR_NONE);
}
