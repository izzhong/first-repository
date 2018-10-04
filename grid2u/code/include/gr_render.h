//file anno
//just a render tool
//this file will not in the engine

#include"gr_clocker.h"
#include<random>
#include<vector>
#include<graphics.h>

#ifndef _GRID_RENDER_H
#define _GRID_RENDER_H

namespace grid
{
	//
	class GridRender
	{
	public:

		//保存窗口信息
		int width_;
		int height_;
		Real fps_;

		//控制渲染帧速
		//这样方便
		Ticker ticker_;
		Timer timer_;
		Ticker test_ticker_;

		char* test_fps;

		std::default_random_engine e;
		std::uniform_int_distribution<int> ui;

	public:

		

	public:

		//initgraph
		GridRender(int width, int height, Real fps, int window_flag);

		//渲染
		void render();

		~GridRender();

	}; //class GridRender

} //namespace grid

#endif //_GRID_RENDER_H