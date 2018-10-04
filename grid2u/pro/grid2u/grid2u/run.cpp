//这里就将是grid2u诞生的地方!

//约束管线已经写好了
//哪怕目前没有碰撞管线应该也可以做一些事情
//所以我决定先试试

#include"game.h"
#include<gr_render.h>

using namespace std;
using namespace grid;
using namespace game;

int main()
{
	//为我的小猫咪创建一个场景
	auto cat = make_shared<LittleCat>();
	auto light = make_shared<Light>();

	//为资源创建资源表
	//如果资源表不是智能指针 仅仅是绑定到这个表上呢??
	auto omt = make_shared<ObjectModelTable>();
	omt->push_front(cat); //这里虽然报错但还是通过了编译 有点奇怪
	omt->push_front(light);

	//资源已经创建好了 为资源创建物理管线
	PhysicsPipeline physics(omt);
	physics.Process();

	//使用图形库创建场景
	GridRender render(1067, 678, 60,SHOWCONSOLE);
	while (true)
	{
		render.render();
	}

	closegraph();
	system("pause");
	return 0;
}