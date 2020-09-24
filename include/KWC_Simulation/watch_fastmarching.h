/*
 * wathch_.fastmarching.h
 *
 *  Created on: Feb 6, 2020
 *      Author: jaekwangkim
 */

#ifndef LOCAL_WATHCH__FASTMARCHING_H_
#define LOCAL_WATHCH__FASTMARCHING_H_


//2d
class WatchFastMarching
{
	public:
  WatchFastMarching(const unsigned int n1,
             const unsigned int n2,
					   const unsigned int lcount,
					   int *labels,
					   double *grid);

	// I am going to connect variables through pointer
	//declared from a constructor
	int *labels; //This will be keep updated from thresholding stage

	FILE *ffmpeg_pipeout;

	//void run ();
	void videoInitialize();
	void videoUpdate_pixel(int updated_index);
	void videoClose();

	private:

  //Video
	int lcount;
	int pcount;
	unsigned char *pixels;

	unsigned int n1;
	unsigned int n2;

	//assert, lcount <256;
	unsigned char *colors;
	double *fast_marching_grid;

};

//Constructor
WatchFastMarching::WatchFastMarching (const unsigned int n1,
                        const unsigned int n2,
						const unsigned int lcount,
						int *labels,
						double *grid)
:n1(n1), n2(n2), lcount(lcount), labels(labels),fast_marching_grid(grid)
{
  pcount= n1*n2;
}

void WatchFastMarching::videoInitialize ()
{
  colors = new unsigned char[lcount+1];
  pixels = new unsigned char[pcount];
  
	assert (lcount<255);

	unsigned int color_gradient = 255 / (lcount+1);

	for(int i=0; i<lcount+1; i++)
	{
		colors[i]= '0'+ (255-color_gradient * i) ;
	}
	std::cout << std::endl;
	
	
	char *string;
	asprintf(&string,"ffmpeg -y -f rawvideo -vcodec rawvideo -pix_fmt gray -s %dx%d -r "
				"60 -i - -f mp4 -q:v 5 -an -vcodec mpeg4 FMout.mp4", n1,n2);

	ffmpeg_pipeout = popen(string, "w");

	//Initial Condition of fast marching
	for(int j=0;j<n2;j++){
	   for(int k=0;k<n1;k++){
		   if(fast_marching_grid[j*n1+k]==0)
		   {pixels[j*n1+k]=colors[labels[j*n1+k]+1];}
		   else
	       {pixels[j*n1+k]=0;}
	   }
    }

	fwrite(pixels, 1, pcount, ffmpeg_pipeout);

}


void WatchFastMarching::videoUpdate_pixel (int updated_index)
{
	pixels[updated_index]=colors[labels[updated_index]+1];

	fwrite(pixels, 1, pcount, ffmpeg_pipeout);
}

void WatchFastMarching::videoClose ()
{
	fflush(ffmpeg_pipeout);
	pclose(ffmpeg_pipeout);

	delete[] pixels;
	delete[] fast_marching_grid;
}




#endif /* LOCAL_WATHCH__FASTMARCHING_H_ */
