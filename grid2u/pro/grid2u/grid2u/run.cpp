//����ͽ���grid2u�����ĵط�!

//Լ�������Ѿ�д����
//����Ŀǰû����ײ����Ӧ��Ҳ������һЩ����
//�����Ҿ���������

#include"game.h"
#include<gr_render.h>

using namespace std;
using namespace grid;
using namespace game;

int main()
{
	//Ϊ�ҵ�Сè�䴴��һ������
	auto cat = make_shared<LittleCat>();
	auto light = make_shared<Light>();

	//Ϊ��Դ������Դ��
	//�����Դ��������ָ�� �����ǰ󶨵����������??
	auto omt = make_shared<ObjectModelTable>();
	omt->push_front(cat); //������Ȼ��������ͨ���˱��� �е����
	omt->push_front(light);

	//��Դ�Ѿ��������� Ϊ��Դ�����������
	PhysicsPipeline physics(omt);
	physics.Process();

	//ʹ��ͼ�οⴴ������
	GridRender render(1067, 678, 60,SHOWCONSOLE);
	while (true)
	{
		render.render();
	}

	closegraph();
	system("pause");
	return 0;
}