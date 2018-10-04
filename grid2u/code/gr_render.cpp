//

#include"include/gr_render.h"
#include<string>

namespace grid
{
	GridRender::GridRender(
		int width, int height, Real fps, int window_flag
	) :
		width_(width), height_(height), fps_(fps),test_fps(new char[6]),
		ticker_(1.0 / fps), timer_(),test_ticker_(0.5),ui(0,0xffffff)
	{
		*test_fps = '\0';
		initgraph(width_, height_, window_flag);
		ticker_.start();
		test_ticker_.start();
		timer_.start();

	}

	void GridRender::render()
	{
		settextcolor(BLUE);
		if (ticker_.tick())
		{					
			BeginBatchDraw();	
			cleardevice();
			
			Real fps = 1.0 / timer_.total_time();
			timer_.restart();		
			if (test_ticker_.tick())
			{
				std::string test_str = std::to_string(fps);
				for (int i = 0; i < 5;++i)
				{
					test_fps[i] = test_str[i];
				}
				test_fps[5] = '\0';
			}
			outtextxy(0, 0, test_fps);			
			EndBatchDraw();
		}
	}

	GridRender::~GridRender()
	{
		delete[] test_fps;
		closegraph();
	}

} //namespace grid