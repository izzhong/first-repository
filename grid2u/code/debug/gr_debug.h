//@header : gr_debug.h
//@abstart : define some debug help function
//@author : zhong
//@date : 2018/8
//@version : bate0.1

#ifndef _GRID_DEBUG_H_
#define _GRID_DEBUG_H_

#include"gr_debug_cout.h"
#include"gr_debug_draw.h"
#include"gr_debug_random.h"

using namespace std;

namespace grid
{
	//@abstract : test gr_math.h
	void TestFileGrMathH()
	{
		//1.21+13.69

		//@test abstract :
		//@date 2018/8/28
		//@test begin
		std::cout << "Test gr_math.h" << std::endl;
		const int size = 1000;
		Real real_array[size];
		Vector2 vector2_array[size];
		Vector3 vector3_array[size];
		RandomRealArray(size, real_array,1,10);
		RandomVector2Array(size, vector2_array,1,10);
		RandomVector3Array(size, vector3_array,1,10);
		CoutRealArray(std::cout, size, real_array);
		CoutVectorArray(std::cout, size, vector2_array);
		CoutVectorArray(std::cout, size, vector3_array);

		Real real_expection = GrExpectation(size, real_array);
		Vector2 vector2_expection = GrExpectation(size, vector2_array);
		Vector3 vector3_expection = GrExpectation(size, vector3_array);
		Cout(std::cout, real_expection) << std::endl;
		Cout(std::cout, vector2_expection) << std::endl;
		Cout(std::cout, vector3_expection) << std::endl;
		Real real_variance = GrVariance(size, real_array);
		Vector2 vector2_variance = GrVariance(size, vector2_array);
		Vector3 vector3_variance = GrVariance(size, vector3_array);
		Cout(std::cout, real_variance) << std::endl;
		Cout(std::cout, vector2_variance) << std::endl;
		Cout(std::cout, vector3_variance) << std::endl;
		//@test end
	}

	void Test()
	{
		//Vector2 rectangle[4];
		//rectangle[0].Set(0.0, 0.0);
		//rectangle[1].Set(100.0, 0.0);
		//rectangle[2].Set(100.0, 100.0);
		//rectangle[3].Set(0.0, 100.0);
		//CoutVectorArray(cout, 4, rectangle);

		//Vector2 unit(1.0, 1.0);
		//unit.Normalize();
		//cout << unit << endl;

		//Matrix3 scale;
		//cout << scale;
		//ScaleMatrix(Vector2(0.0, 0.0), unit, -1.0, &scale);
		//cout << scale;

		//scale.Transform(4, rectangle);
		//CoutVectorArray(cout, 4, rectangle);

		//cout << GrCos(45) << std::endl;
		//cout << GrCos(MATH_CONSTANT_PI / 4.0) << endl;
		//cout << GrCos(AngleToRadian(45)) << endl;

		//Matrix3 rotate;
		//RotateMatrix(Vector2(0.5, 0.5), 45.0, &rotate);
		//rotate.Transform(4, rectangle);
		//CoutVectorArray(cout, 4, rectangle);

		//initgraph(1067, 600);
		//setorigin(534, 300);
		//setaspectratio(1.0, -1.0);

		//POINT rectangle_point[4];
		//for (int i = 0; i < 4; ++i)
		//{
		//	rectangle_point[i].x = static_cast<long>(rectangle[i].x);
		//	rectangle_point[i].y = static_cast<long>(rectangle[i].y);
		//}
		//
		//while (1)
		//{
		//	BeginBatchDraw();
		//	line(0, -1000, 0, 1000);
		//	line(-1000, 0, 1000, 0);
		//	polygon(rectangle_point, 4);
		//	EndBatchDraw();
		//}
	}

	//@date : 2018/9/1
	//@result	:	大多数函数在一般情况下工作良好
	//				但是存在浮点数精度问题
	//				当浮点数非常小时 计算就会出错
	//				后面测试矩阵变换时这个问题可能更加严重
	void TestFileGrVectorH()
	{
		Vector2 vec21(61,-3.6);
		Vector2 vec22(1.0);
		Vector2 vec23(2.0, 3.5);
		Vector3 vec31(-21,5.1,5);
		Vector3 vec32(-2);
		Vector3 vec33(31, 31, 3);

		Cout(std::cout, vec21) << std::endl;
		Cout(std::cout, vec22) << std::endl;
		Cout(std::cout, vec23) << std::endl;
		Cout(std::cout, vec31) << std::endl;
		Cout(std::cout, vec32) << std::endl;
		Cout(std::cout, vec33) << std::endl;
		Cout(std::cout, vec21 + vec22);
		Cout(std::cout, vec21 - vec22);
		Cout(std::cout, vec21 * vec22);
		Cout(std::cout, vec21 += vec22);
		Cout(std::cout, vec21 -= vec22);
		Cout(std::cout, vec21 *= vec22);
		Cout(std::cout, vec21 * MATH_CONSTANT_PI);
		Cout(std::cout, vec21 *= MATH_CONSTANT_E);
		Cout(std::cout, vec21 / 3.14);
		Cout(std::cout, vec21 /= -3.15);
		Cout(std::cout, vec21 == vec22);
		Cout(std::cout, vec21 != vec22);

		Cout(std::cout, vec31 + vec33);
		Cout(std::cout, vec31 - vec33);
		Cout(std::cout, vec31 * vec33);
		Cout(std::cout, vec31 += vec33);
		Cout(std::cout, vec31 -= vec33);
		Cout(std::cout, vec31 *= vec33);
		Cout(std::cout, vec31 * MATH_CONSTANT_PI);
		Cout(std::cout, vec31 *= MATH_CONSTANT_E);
		Cout(std::cout, vec31 / 3.14);
		Cout(std::cout, vec31 /= -3.15);
		Cout(std::cout, vec31 == vec33);
		Cout(std::cout, vec31 != vec33) << std::endl;

		Cout(std::cout, "------------------------------------------------------------------") << std::endl;

		Cout(std::cout, vec21);
		Cout(std::cout, vec22);
		Cout(std::cout, vec23);
		Cout(std::cout, vec21.DotProduct(vec22)) << std::endl;;
		Cout(std::cout, vec22.CrossProduct(vec23)) << std::endl;

		Cout(std::cout, vec31);
		Cout(std::cout, vec32);
		Cout(std::cout, vec33);
		Cout(std::cout, vec31.DotProduct(vec33)) << std::endl;
		Cout(std::cout, vec32.CrossProduct(vec33));

		Cout(std::cout, "------------------------------------------------------------------") << std::endl;
		//cov
		const int size = 3;
		Vector2 vector2_array[size]{};
		RandomVector2Array(size, vector2_array, 0, 600);
		CoutVectorArray(std::cout, size, vector2_array);
		std::cout << "expection : " << GrExpectation(size, vector2_array) << std::endl;
		std::cout << "covariance : " << Covariance(size, vector2_array) << std::endl;
		std::cout << "cov coe : " << CorrelationCoefficient(size, vector2_array) << std::endl;

		Vector2 A(0.0, 0.0);
		Vector2 B(1.0, 0.0);
		Vector2 C(0.0, 1.0);
		Vector2 AB = B - A;
		Vector2 AC = C - A;
		Vector2 triple = VectorTripleProduct(AB, AC, AB);
		Vector2 invert_triple;
		Real radian = 0.0;
		Real angle = 0.0;

		//VectorAngle单独测试
		A.Set(1.0, 0.0);
		B.Set(0.0, 1.0);
		cout << A << B;
		radian = VectorRadian(A, B);
		cout << "--------------------------------------------" << endl;
		std::cout << radian << std::endl;
		std::cout << RadianToAngle(radian) << std::endl;

		A.Set(0.000001, 0.000001);
		B.Set(0.0, 1.0);
		cout << A << B;
		radian = VectorRadian(A, B);
		cout << "--------------------------------------------" << endl;
		std::cout << radian << std::endl;
		std::cout << RadianToAngle(radian) << std::endl;

		A.Set(1.0, 0.0);
		B.Set(-1.0, 0.0);
		cout << A << B;
		radian = VectorRadian(A, B);
		cout << "--------------------------------------------" << endl;
		std::cout << radian << std::endl;
		std::cout << RadianToAngle(radian) << std::endl;

		A.Set(0.0, 0.0);
		B.Set(0.0, 1.0);
		cout << A << B;
		radian = VectorRadian(A, B);
		cout << "--------------------------------------------" << endl;
		std::cout << radian << std::endl;
		std::cout << RadianToAngle(radian) << std::endl;

		A.Set(1.0, 0.0);
		B.Set(-2.0, 1.0);
		cout << A << B;
		radian = VectorRadian(A, B);
		cout << "--------------------------------------------" << endl;
		std::cout << radian << std::endl;
		std::cout << RadianToAngle(radian) << std::endl;

		//我来指定一个圆 让这个函数来无死角的计算各个角度
		Real y = 0.0;
		Vector2 OA(1.0, 0.0);
		Real r = 0.0001;
		cout << "circle test start" << endl;
		for (Real x = -r; x <= r; x += 0.1 * r)
		{
			y = GrSqrt(1 - x * x);
			B.Set(x, y);
			C.Set(x, -y);
			cout << B << C;
			radian = VectorRadian(OA, C);
			cout << "--------------------------------------------" << endl;
			std::cout << radian << std::endl;
			std::cout << RadianToAngle(radian) << std::endl;
		}

		//相关系数测试
		default_random_engine e(static_cast<unsigned int>(time(0)));
		uniform_real_distribution<Real> ur(-1.0, 1.0);
		Vector2 cov[100];
		for (int i = 0; i < 100; ++i)
		{
			cov[i].Set(i + ur(e), i + ur(e));
		}
		std::cout << CorrelationCoefficient(100, cov) << endl;

		//初步测试成功
		//下一步进行正确性和健壮性测试
		//首先实现向量可视化
		initgraph(1067, 600, SHOWCONSOLE);
		setorigin(1067 / 2, 600 / 2);

		while (true) 
		{
			BeginBatchDraw();
			cleardevice();		
			setfillcolor(0x49CFA7);
			fillrectangle(-1000, 1000, 1000, -1000);
			setlinecolor(BLACK);
			line(0, -1000, 0, 1000);
			line(-1000, 0, 1000, 0);
			setlinecolor(BLUE);
			RandomVector2Array(size, vector2_array, -500.0, 500.0, -300, 300);
			for (int i = 0; i < size; ++i)
			{
				DrawVector2(vector2_array[i]);
				DrawString(vector2_array[i], std::to_string(i + 1));
			}
			A = vector2_array[0];
			B = vector2_array[1];
			C = vector2_array[2];
			AB = B - A;
			AC = C - A;	
			triple = VectorTripleProduct(AB, AC, AB);
			invert_triple = VectorTripleProduct(AB, -AC, AB);
			radian = VectorRadian(triple, AB);
			angle = RadianToAngle(radian);
			std::cout << "----------------------------------" << std::endl;
			std::cout << "AB : " << AB;
			std::cout << "tri : " << triple;
			std::cout <<"tri * AB : "<< triple.DotProduct(AB) << std::endl;
			std::cout << "tri _|_ AB : " << IsVectorVertical(AB, triple) << std::endl;
			std::cout << "angle(tri,AB) : " << angle << std::endl;
			std::cout << "tri x inv_tri :" << triple.CrossProduct(invert_triple) << std::endl;
			std::cout << "tri || inv_tri :" << IsVectorParallel(triple, invert_triple) << std::endl;
			std::cout << "angle(1,2) : " << RadianToAngle(VectorRadian(A, B)) << std::endl;
			std::cout << "angle(2,3) : " << RadianToAngle(VectorRadian(B, C)) << std::endl;
			std::cout << "angle(3,1) : " << RadianToAngle(VectorRadian(C, A)) << std::endl;
			DrawVector2(A,B);
			if (!triple.IsZero()) {
				triple.Normalize();
			}
			DrawVector2((A+B)*0.5, (A + B)*0.5 + triple * 100);
			DrawString((A + B)*0.5 + triple * 100, "tri");
			Sleep(1000);
			EndBatchDraw();		
		}


	}

	//@date : 2018/9/2
	//@result : 初步测试工作量好
	//			未测试变换逆矩阵
	//			未测试协方差矩阵
	void TestFileGrMatrixH()
	{
		//常规测试
		Matrix2 m2;
		Matrix3 m3;
		m2.Set(1.0, 1.0, 1.0, 1.0);
		m3.Set(1, 1, 1, 1, 1, 1, 1, 1, 1);
		m2 <<= 1, 2, 3, 4;
		m3 <<= 1, 2, 3, 4, 5, 6, 7, 8, 9;
		m2 <<= 5, 6, 7, 8;
		m3 <<= 9, 8, 7, 6, 5, 4, 3, 2, 1;
		TestMatrixCout(m2 + m2);
		TestMatrixCout(m3 - m3);

		//部分函数的正确性测试
		m2 <<= 3, -2, 2, 1;
		cout << m2.Determinant() << endl;
		m3 <<= 1, 2, -4, -2, 2, 1, -3, 4, -2;
		TestMatrixCout(m3);
		cout << m3.Determinant() << endl;
		Matrix2 A(-2, 4, 1, -2);
		Matrix2 B(2, 4, -3, -6);
		cout << A * B;
		cout << B * A;
		Matrix3 A3(1, 1, 1, 1, 1, -1, 1, -1, 1);
		Matrix3 B3(1, 2, 3, -1, -2, 4, 0, 5, 1);
		cout << A3 * B3 * 3 - A3 * 2;
		cout << A3.Transpose() * B3;
		m3 <<= 1, 2, 3, 2, 2, 1, 3, 4, 3;
		cout << m3.Invert();
		A3.Set(4, 3, 1, 1, -2, 3, 5, 7, 0);
		Vector3 v3(7, 2, 1);
		cout << A3 * v3 << endl;
		A.Set(1, 2, 2, 5);
		A3 <<= 1, 2, -1, 3, 4, -2, 5, -4, 1;
		cout << A.Invert() << endl;
		cout << A3.Invert() << endl;
		Matrix3 E;
		A3 <<= 0, 3, 3, 1, 1, 0, -1, 2, 3;
		cout << (A3 - E * 2).Invert() * A3;
		A <<= 2, 5, 1, 3;
		B <<= 4, -6, 2, 1;
		cout << A.Invert() * B;
		Vector2 v2(-2, 3);
		A <<= 1, 2, 3, 4;
		cout << A * v2 << endl;

		//特征值与特征向量
		m2 <<= 3, -1, -1, 3;
		Vector2 ev1, ev2;
		m2.EnginVector(&ev1, &ev2);
		Real er1, er2;
		m2.EnginValue(&er1, &er2);
		cout << ev1 << ev2;
		cout << er1 << "," << er2 << endl;

		//矩阵变换测试
		Vector2 translate_vector(100,100);
		Matrix3 translate_matrix;
		TranslateMatrix(translate_vector, &translate_matrix);

		constexpr int size = 100;
		Vector2 vectors[size];

		//协方差矩阵


		initgraph(1067, 600, SHOWCONSOLE);
		setorigin(1067 / 2, 600 / 2);

		while (true)
		{
			BeginBatchDraw();
			cleardevice();
			setfillcolor(0x49CFA7);
			fillrectangle(-1000, 1000, 1000, -1000);
			setlinecolor(WHITE);
			line(0, -1000, 0, 1000);
			line(-1000, 0, 1000, 0);
			setlinecolor(BLUE);
			RandomVector2Array(size, vectors, -500.0, 500.0, -300, 300);
			for (int i = 0; i < size; ++i)
			{
				DrawVector2(vectors[i]);
				DrawString(vectors[i], std::to_string(i + 1));
			}

			TranslateTransform(
				-vectors[0] * 2, size, vectors
			);
			setlinecolor(RED);
			for (int i = 0; i < size; ++i)
			{
				DrawVector2(vectors[i]);
				DrawString(vectors[i], std::to_string(i + 1));
			}		

			RotateTransform(
				Vector2(0.0, 0.0), -72, size, vectors
			);
			setlinecolor(YELLOW);
			for (int i = 0; i < size; ++i)
			{
				DrawVector2(vectors[i]);
				DrawString(vectors[i], std::to_string(i + 1));
			}

			ProjectTransform(Vector2(0, 0), Vector2(0, 1), size, vectors);
			setlinecolor(0x262626);
			for (int i = 0; i < size; ++i)
			{
				DrawVector2(vectors[i]);
				DrawString(vectors[i], std::to_string(i + 1));
			}

			ReflectTransform(Vector2(0, 0), Vector2(0, 1), size, vectors);
			setlinecolor(BROWN);
			for (int i = 0; i < size; ++i)
			{
				DrawVector2(vectors[i]);
				DrawString(vectors[i], std::to_string(i + 1));
			}
			EndBatchDraw();
			Sleep(1000);
		}
	}

	//@date : 2018/9/2
	void TestFileGrGeometryH()
	{
		//point set
		const int size = 5;
		Vector2 vec[size];
		RandomVector2Array(size, vec);
		PointSet ps(size,vec);
		cout << ps;
		for (int i = 0; i < size; ++i)
		{
			ps.EraseHead();
			cout << ps;
		}
		RandomVector2Array(size, vec);
		ps.Set(size, vec);
		cout << ps;
		for (int i = 0; i < size; ++i)
		{
			ps.EraseTail();
			cout << ps;
		}
		ps.Set(size, vec);
		cout << ps;
		//ps.Clear();
		cout << ps;
		
		auto p = ps.Find(
			[](const Point& point) {
			return point.x > 0.5; 
		});

		cout << "[erase ptr]" << p << endl;
		ps.Erase(p);
		cout << ps;

		p = ps.Find(
			[](const Point& point) {
			return point.x > 0.5;
		});
		auto ip = ps.InsertAfter(p, Vector2(233, 233));
		cout << " [insert ptr]" << p << endl;
		cout << " [inserr after ptr] " << ip << endl;
		cout << ps;

		p = ps.Find(
			[](const Point& point) {
			return point.x < 0.5;
		});
		ip = ps.InsertBefore(p, Vector2(233, 233));
		cout << " [insert ptr]" << p << endl;
		cout << " [inserr before ptr] " << ip << endl;
		cout << ps;
		ps.ForEach(CoutPoint);

		cout << ps;
		Matrix3 translate;
		TranslateMatrix(Vector2(100, 100), &translate);
		ps.Transform(translate);
		cout << ps;

		//@date : 2018/9/8
		//bug merge与std::move一起工作会出现异常 
		//		目前我们仍不知道这个bug的原因是什么
		cout << "--------------------------------------" << endl;
		PointSet s1;
		int size1 = 10;
		PointSet s2;
		int size2 = 5;
		for (int i = 0; i < size1; ++i)
		{
			s1.Push(Point(i, i));
		}
		for (int i = 0; i < size2; ++i)
		{
			s2.Push(Point(i, i));
		}
		cout << s1;
		cout << s2;
		s1.Merge(s2);
		cout << s1;
		cout << s2;
		PointSet&& s3 = std::move(s2);
		//s1 = std::move(s2);
		//s1.Merge(std::move(s2));
		s1.Merge(s3);
		s1.Clear();
		s2.Clear();
		s3.Clear();
		cout << s1;
		cout << s2;
		cout << s3;

		//PointArray
		//仍然有这个bug与std::move协同工作会出现问题
		PointArray pa1{ {1,1},{2,2},{3,3} };
		PointArray pa2(pa1);
		CoutVectorArray(cout, pa1.Size(), pa1.Array());
		CoutVectorArray(cout, pa2.Size(), pa2.Array());
		pa1 = pa2;
		CoutVectorArray(cout, pa1.Size(), pa1.Array());
		CoutVectorArray(cout, pa2.Size(), pa2.Array());
		for (int i = 0; i < pa1.Size(); ++i)
		{
			cout << pa1[i];
		}


		//Plygon Test

		Polygon poly1{ {0,0},{100,0} ,{100,100},{0,100} };
		Polygon polycir(Circle({0,0}, 100),5);
		Polygon polyarc(Arc({ 0,0 }, { 200,0 }, { 0,200 }, 200),3);
		Polygon polysec(Sector({ 0,0 }, { -200,0 }, { 0,200 }, 200),3);
		Polygon polyseg(Segment({ 0,0 }, { 300,300 }));
		Polygon polyray(Ray({ 0,0 }, { -0.707,-0.707 }), 10000);
		Polygon polyrec(Rectangle(100, -100, 150, 150));
		Polygon polytri(Triangle({ {0,0},{-100,0},{-100,-100} }));
		Polygon polybox(Box({ 0,0 }, { 0.707,0.707 }, { -0.707,0.707 }, 300, 100));
		initgraph(1067, 600);
		setorigin(1067 / 2, 600 / 2);
		Real angle = 0.0;
		Matrix3 rotate;
		while (true)
		{
			BeginBatchDraw();
			cleardevice();
			//axis of x-y
			setlinecolor(GREEN);
			line(0, -1000, 0, 1000);
			line(1000, 0, -1000, 0);
			setlinecolor(WHITE);
			RotateMatrix(Point(0, 0), angle+=0.001, &rotate);
			poly1.Transform(rotate);
			polycir.Transform(rotate);
			polyarc.Transform(rotate);
			polysec.Transform(rotate);
			polyseg.Transform(rotate);
			polyray.Transform(rotate);
			polybox.Transform(rotate);
			DrawPolygon(poly1);
			DrawPolygon(polycir);
			DrawPolygon(polyarc);
			DrawPolygon(polysec);
			DrawPolygon(polyseg);
			DrawPolygon(polyray);
			DrawPolygon(polyrec);
			DrawPolygon(polytri);
			DrawPolygon(polybox);
			EndBatchDraw();
		}
		closegraph();
	}

	//@date : 2018/9/8
	void TestFileBoundingH()
	{
		int size = 1000;
		PointArray pa(size);	
		PointSet ps;
		AABB aabb;
		CBC cbc;
		OBB obb;
		Circle cir;
		Box box;
		Triangle tri;
		
		initgraph(1067, 600,SHOWCONSOLE);
		setorigin(1067 / 2, 600 / 2);
		while (true)
		{
			BeginBatchDraw();
			cleardevice();
			setlinecolor(WHITE);
			RandomVector2Array(
				size, pa.Array(), -300, 300, -200, 200
			);
			for (int i = 0; i < pa.Size(); ++i)
			{
				DrawPoint(pa[i]);
			}

			obb.CreateBounding(pa.Size(), pa.Array());
			setlinecolor(WHITE);
			DrawBox(obb.GetGeometry());
			setlinecolor(RED);
			DrawVector2(obb.GetGeometry().center, obb.GetGeometry().center + obb.GetGeometry().direction_x * 100);
			setlinecolor(GREEN);
			DrawVector2(obb.GetGeometry().center, obb.GetGeometry().center + obb.GetGeometry().direction_y * 100);

			EndBatchDraw();
		}
	}

	//@date : 2018/9/22
	void TestFileIntersectionH()
	{
		Point p1;
		Point p2(1,1);
		cout << IntersectionTest(p1, p2) << endl;

	}

} //namespace grid

#endif //_GRID_DEBUG_H_