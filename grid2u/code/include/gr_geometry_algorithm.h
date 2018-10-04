//@header : gr_grometry_algorithm.h
//@author : zhong
//@date	  : 2018/9

//TODO(zhong) : 相交测试写成相同的样子有助于管线的模板化
//				bool IntersectionTest(const Geometry left,const Geometry& right);

#ifndef _GRID_GEOMETRY_ALGORITHM_H_
#define _GRID_GEOMETRY_ALGORITHM_H_

#include"gr_geometry.h"

namespace grid
{
	
	//------直线线段射线算法----------

	void TwoPointsLineEquation(
		const Point& A, const Point& B,
		Real* coeA, Real* coeB, Real* coeC
	);

	void PointNormalLineEquation(
		const Point& P,
		const Vector2& normal,
		Real* coeA, Real* coeB, Real* coeC
	);


	//其实这个是线段跨立试验
	bool SegmentStraddeleTest(
		const Segment& segment_1,
		const Segment& segment_2
	);

	////若A B在直线CD两侧 返回真
	bool LineStraddleTest(
		const Point& A,
		const Point& B,
		const Point& LineC,
		const Point& LineD
	);

	bool LineStraddleTest(
		const Point& A,
		const Point& B,
		const Line& line
	);

	//判断向量 OC 是否在向量 OA OB 之间
	bool VectorRegionTest(
		const Vector2& OA, 
		const Vector2& OB,
		const Vector2& OC
	);

	//-----最远点算法------

	////获得在这个方向上最远点的下标
	int GetFartestPointIndex(
		const Vector2 direction,
		int size, const Point* point_array
	);

	Point GetFartestPoint(
		const Vector2 direction,
		int size, const Point* point_array
	);

	PointNode* GetFartestPointPtr(
		const Vector2& direction,
		const PointSet& pointset
	);

	Point GetFartestPoint(
		const Vector2& direction,
		const PointSet& pointset
	);

	//-----最近点------

	Point ClosestPointPtSeg(
		const Point& point, const Segment& segment
	);

	Point ClosestPointPtRec(
		const Point& point,const Rectangle& rectangle
	);

	Point ClosestPointPtTri(
		const Point& point, const Triangle& triangle
	);

	Point ClosestPointPtBox(
		const Point& point, const Box& box
	);

	Point ClosestPointPtPol(
		const Point& point, const Polygon& polygon
	);

	Point ClosestPointSegSeg(
		const Segment& segment_1, const Segment& segment_2
	);

	Point ClosestPointRayRay(
		const Ray& ray_1, const Ray& ray_2
	);

	Point ClosestPointSegTri(
		const Segment& segment, const Triangle& triangle
	);

	Point closestPointTriTri(
		const Triangle& triangle_1, const Triangle& triangle_2
	);

	//-----距离算法-------

	////@abstract : calculate the distance between two points
	Real CalculatePointsDistance(
		const Point& point_1,
		const Point& point_2
	);

	//@abstract : calculate the cloest distance between 
	//			  a point and a segment
	Real CalculatePointToSegmentDistance(
		const Point& point,
		const Segment& segment
	);

	//----投影计算算法-------
	//@abstract : 计算投影长度
	Real CalculateProjectionLength(
		const Vector2& axis,
		int size, const Point* point_array,
		Real* min, Real* max,
		int* min_nptr = nullptr, int* max_nptr = nullptr
	);

	//传出参数可以给空指针 空指针就不做处理就可以了
	Real CalculateProjectionLength(
		const Vector2& axis,
		const PointSet& pointset,
		Real* min, Real *max,
		PointNode* min_nptr = nullptr, PointNode* max_nptr = nullptr
	);

	Real CalculateProjectionLength(
		const Vector2& axis,
		const Circle& circle,
		Real* min, Real* max
	);

	Real CalculateProjectionLength(
		const Vector2& axis, const Arc& arc,
		Real* min, Real* max
	);

	Real CalculateProjectionLength(
		const Vector2& axis, const Sector& sector,
		Real* min, Real* max
	);

	Real CalculateProjectionLength(
		const Vector2& axis, const Triangle& triangle,
		Real* min, Real* max
	);

	Real CalculateProjectionLength(
		const Vector2& axis, const Rectangle& rectangle,
		Real* min, Real* max
	);

	Real CalculateProjectionLength(
		const Vector2& axis, const Box& box,
		Real* min, Real* max
	);

	Real CalculateProjectionLength(
		const Vector2& axis, const Polygon& polygon,
		Real* min, Real*max
	);

	//包围体算法
	Rectangle RectangleBounding(const Polygon& polygon);

	//--------多边形相关算法-------

	void FindHull(
		const PointSet& Sk,
		PointNode* P, PointNode* Q,
		PointSet* polygon
	);

	void QuickHull(
		int size, Point* point_array,
		PointSet* pointset
	);

	////实现凸多边形检测算法
	///*
	//非常简单
	//凸多边形就是以任意一条边画一条线 剩下的所有点都在这条线的一端
	//所以判断点与线的方向关系 用点乘 或者 叉乘都可以
	//叉积的一个非常重要性质是可以通过它的符号判断两矢量相互之间的顺逆时针关系
	//算法
	//顺序选取三个点 这三个点可以确定这个多边形的旋转方向
	//剩下的依次向下选取点 一旦旋转方向与初始方向不同
	//那他就是一个凹多边形
	//否则就是凸多边形
	//算法时间复杂度 o[n]
	//*/
	bool IsConvexPolygon(
		const Polygon& polygon
	);

	//------------------------------------------------


	//相离检测算法

	template<typename Geometry1,typename Geometry2>
	bool SeparationDetection(
		const Geometry1& gry_1, const Geometry& gry_2
	);

	bool CircleSeparationDetection(
		const Circle& circle_1,
		const Circle& circle_2
	);

	bool RectangleSeparationDetection(
		const Rectangle& rectangle_1,
		const Rectangle& rectangle_2
	);

	bool BoxSeparationDetection(
		const Box& box_1,
		const Box& box_2
	);

	bool PolygonSeparationDetection(
		const Polygon& polygon_1,
		const Polygon& polygon_2
	);

	//相交检测

	bool IntersectionDetection();

	//------碰撞检测算法----

	//几何相交状态
	enum class GeometryIntersectionState : __int8
	{
		SEPERATION,
		OUTER_CONTACT,
		INTERSECTION,
		INTER_CONTACT,
		CONTAINT
	};

	struct CollisionState
	{
		Vector2 collision_vector;
		Real collision_depth;
		Point closest_point[2];
	};

	//几何学相交检测算法结果
	//没有必要每个碰撞都检测如此复杂
	//大部分时候我们只是想直到他们有没有碰撞
	struct IntersectionResult
	{
		GeometryIntersectionState state;
		Vector2	intersection_vector;
		Real	intersection_depth;
		Point	closet_point[2];//最近点就是碰撞点		
	};

	bool IsPointSegmentCollineation(
		const Point& point,
		const Segment& segment
	);

	bool IsPointOnSegment(
		const Point& point,
		const Segment& segment
	);

	bool IsPointInTriangle(
		const Point& point,
		const Triangle& triangle
	);

	bool IsPointOnTriangle(
		const Point& point,
		const Triangle& triangle
	);

	bool IsPointOutTriangle(
		const Point& point,
		const Triangle& triangle
	);

	////rectangle
	///*
	//不可旋转的矩形的碰撞检测
	bool IsRectangleIntersection(
		const Rectangle& rectangle_1,
		const Rectangle& rectangle_2
	);

	//point to point intersection detection
	void PointToPointIntersectionDetection(
		const Point& point_1,
		const Point& point_2,
		IntersectionResult *result
	);

	//对简单几何图形直接应用EPA算法 
	//计算出公式直接计算就可以了 
	//比这些向量计算什么的快多了
	//计算点与线段的最短距离 和 最短距离点
	void PointToSegmentIntersectionDetection(
		const Point& point,
		const Segment& segment,
		IntersectionResult* result
	);

	void PointToTriangleIntersectionDetection(
		const Point& point,
		const Triangle& triangle,
		IntersectionResult* result
	);

	void PointToCircleIntersectionDetection(
		const Point& point,
		const Circle& circle,
		IntersectionResult* result
	);
	
	void PointToRectangleIntersectionDetection(
		const Point& point,
		const Rectangle& rectangle,
		IntersectionResult* result
	);

	void PointToBoxIntersectionDetection(
		const Point& point,
		const Box& box,
		IntersectionResult* result
	);

	void PointToSectorIntersectionDetection(
		const Point& point,
		const Sector& sector,
		IntersectionResult* result
	);

	void PointToArcIntersectionDetection(
		const Point& point,
		const Arc& arc,
		IntersectionResult* result
	);

	////判断点是否在多边形内
	//inline bool point_in_polygon(Vector2 _point, polygon _pol)
	//{
	//	//射线法
	//	//从无限远处取一点 向着指定点移动
	//	//如果遇到多边形的一点则算作进入多边形 
	//	//如果再遇到一点 表示出了多边形

	//	//总之遇到奇数个点代表点在多边形内
	//	//遇到偶数个点表示不在多边形内

	//	//其实真正实现起来大概就是求以指定点位一点 任意向量位线向量的直线与多边形所有边的交点把

	//	//针对某些特殊情况 在对多边形的每条边进行判断时 要考虑一下特殊情况
	//	/*
	//	假设当前处理的边是p1p2
	//	1.如果点在p1p2边上 则直接判定点p在多边形内
	//	2.如果从p发出的射线正好穿过p1或p2 那么这个交点会被算作两次
	//	对这种情况的处理原则是: 如果p的y坐标与p1,p2中较小的y坐标相同 则忽略这个交点
	//	3.如果从p发出的射线与p1 p2 平行 则忽略这条边
	//	*/

	//	/*
	//	真正的工程代码 还会增加一个多边形的外包矩形快速判断 对不再多边形周围的情况做快速排出
	//	求多边形的外包矩阵算法
	//	遍历多边形的所有节点找出各个方向上的最值即可
	//	*/
	//}
	void PointToPolygonIntersectionDetection(
		const Point& point,
		const Polygon& polygon,
		IntersectionResult* result
	);

	////判断线段是否相交
	//分两个步骤完成
	//1.快速排斥试验
	//设以线段p1p2为对角线的矩形为R1,设以线段Q1Q2为对角线的矩形为R2
	//如果R1 R2不相交 则两线段不会有交点
	//2.跨立试验
	//如果两线段相交 则两线段必然相互跨立对方
	//所谓跨立 指的是一条线段的两断电分别位于另一线段所在直线的两边
	//若p1p2跨立q1q2 则矢量 p1-q1 和 p2-q2位于矢量q2-q1的两侧
	//所以有
	//(p1-q1)x(q2-q1) * (p2-q1)x(q2-q1) < 0
	void SegmentIntersectionDetection(
		const Segment& segment_1,
		const Segment& segment_2,
		IntersectionResult* result
	);

	void SegmentRectangleIntersectionDetection(
		const Segment& segment,
		const Rectangle rectangle,
		IntersectionResult* result
	);

	//超级提炼总结
	//线段与其他图形的相交检测可以这样计算
	/*
		首先执行快速检测 矩形检测
		然后执行跨立检测
	*/

	////直线检测
	////射线和圆形的检测
	//inline bool line_circle(
	//	Vector2 _line_v, Vector2 _point, circle _cir
	//)
	//{
	//	//判断点到直线的距离与半径的差值

	//	//射线需要判断 是否当直线相交时 但是射线并不相交 
	//	//只需要计算 从发射点到圆心的的向量与直线向量的向量积与0的关系即可
	//}
	void SegmentCircleIntersectionDetection(
		const Segment& segment,
		const Circle& circle,
		IntersectionResult* result
	);

	
	//圆形和矩形
	//inline bool circle_rectangle(
	//	cir* _cir, rectangle _rec
	//)const
	//{
	//	Vector2 close;
	//	//首先确定圆心与矩形上最近的点
	//	if (_cir->x <= _rec->x) close.x = _rec->x;
	//	else if (_cir->x >= _rec->x + _rec->width) close.x = _rec->x + _rec->width;
	//	else close.x = _cir->x;
	//	//y也一样

	//	//然后计算此点到圆心的距离并与半径比较
	//	return (close - Vector2(_cir->x, _cir->y)).lengthSquare() < _cir->radius * _cir->radius;
	//}
	void RectangleCircleIntersectionDetection(
		const Rectangle& rectangle,
		const Circle& circle,
		IntersectionResult* result
	);

	////圆形与旋转矩形
	///*
	//虽然矩形以其中心为轴进行了旋转 但是判断它与圆形是否发生了碰撞本质上还是找出矩形上离圆心的最近点
	//我们可以把整个图形看作对画布进行了依次旋转
	//使之变成常规的判断
	//*/
	void BoxCircleIntersectionDetection(
		const Box& box,
		const Circle& circle,
		IntersectionResult* result
	);

	//矩形与旋转盒
	void RectangleBoxIntersectionDetection(
		const Rectangle& rectangle,
		const Box& box,
		IntersectionResult* result
	);

	//旋转盒碰撞检测
	void BoxIntersectionDetection(
		const Box& box_1,
		const Box& box_2,
		IntersectionResult* result
	);

	void CircleIntersectionDetection(
		const Circle& circle_1,
		const Circle& circle_2,
		IntersectionResult* result
	);



	////添加一组扇形碰撞的算法
	////扇形和圆形的碰撞判断
	//bool sector_circle()
	//{
	//	Vector a(sector.vertex[0] - sector.center);
	//	Vector b(sector.vertex[1] - sector.center);
	//	Vector c(circle.center - sector.center);
	//	//判断圆心是否在扇形区域内
	//	if (a x c * c x b > 0 && a x c * a x b > 0)
	//	{
	//		//直接使用圆形碰撞检测算法
	//	}
	//	else
	//	{
	//		//计算圆心与线段的最短距离
	//		//与半径作比较即可
	//	}
	//}
	//圆形和扇形
	//这个东西的检测十分复杂 并不像上面所述这么简单
	//暂时放弃实现
	void CircleSectorIntersectionDetection(
		const Circle& circle,
		const Sector& sector,
		IntersectionResult* result
	);


	////SAT
	bool SAT(
		const Polygon& polygon_1,
		const Polygon& polygon_2
	);


	Vector2 GetClosetPointOnPolygon(
		const Point& point,
		const Polygon& polygon
	);



	////圆形与多边形的碰撞检测
	///*
	//只需要将圆形投射到一条投影轴上即可
	//这条轴就是圆心与多边形定点中最近的一点的连线
	//*/
	bool CirclePolygonSAT(
		const Circle& circle,
		const Polygon& polygon
	);



		////-------------------------------------------------------------------
		///*
		//非常牛逼的一部分内容
		//用于计算相交深度
		//GJK算法

		//!!:如果两个形状重叠或者相交那么这两个形状的闵科夫斯基差产生的点的外围形状会包含原点
		//*/

		////GJK算法支持函数
		//Vector2 support(polygon _poy1, polygon _poy2, Vector2 _dir)
		//{
		//	//获得两个多边形在给定方向上最远的点
		//	//其实就是找到两个多边形在这个方向上相距最远的两个点
		//	//作为一个名科夫斯基差
		//	//简单 用数量积就可以了

		//	//那么怎么获得这个多边形在这个方向上的最远的点呢??
		//	for (size_t i = 0; i < _poy1.side; ++i)
		//	{
		//		//计算各个顶点在给定向量上的分量 取得最大值
		//		//在计算过程中 保存最大值节点的下标
		//	}
		//	//对于poy2 则使用-dir进行计算 取得最大分量
		//	for (size_t i = 0; i < _poy2.side; ++i)
		//	{
		//		//计算各个顶点在给定向量上的分量 取得最大值
		//		//在计算过程中 保存最大值节点的下标
		//	}

		//	//返回这两个点的名科夫斯基差

		//}

		////迭代方案
		//void GJK()
		//{
		//	vector d = ; //最初的方向选择也对迭代方案产生一些影响
		//				 //当多边形点的数目较少时 我们可以采取一些方案
		//				 //什么点最能代表一个多边形呢 当然是中心!!所以示例代码中使用了以两个多边形中心为起点终点的向量作为初始向量
		//				 //我总觉得应该是垂直向量
		//				 ///如果两个多边形重合 中心就一样 可以直接判定相交
		//				 //我要试一下
		//	simplex.add(support(A, B, d));
		//	d.negate();
		//	while (true)
		//	{
		//		simplex.add(support(A, B, d));
		//		//检测刚刚添加的两个点是否跨过了原点 
		//		if (simplex.getLast().dot(d) <= 0)
		//		{
		//			//判断 根据现在 的 d 生成的名科夫斯基差点
		//			//是否在这个方向上跨过了原点
		//			//我们从原点 与d垂直的方向画一条直线
		//			//如果跨过了原点 那么指向将把这两个的连线分开 
		//			//也就是 点乘 >0
		//			//如果点乘为0 那么新生成的点就在直线上
		//			//如果点乘<0 那么新生成的点就没有在这个方向跨国原点 
		//			//这里还需要研究 == 0到底是什么情况 是接触吗??


		//			return false;
		//		}

		//		else
		//		{
		//			//现在我们需要确定原点是否在单纯行中
		//			if (containsOrigin(simplex, d))
		//			{
		//				//如果发生了碰撞
		//				return true;
		//			}

		//		}
		//	}
		//}

		////在计算时  
		////如果多边形点接触 那么会有一个法向量在计算时会计算成零向量
		////
		////还有很多情况需要分析啊 单纯的边接触  包含  包含点接触  包含边接触 。。。。

		////可以优化 
		////如果我们计算出了一个单纯形的点 就是原点 那么他们肯定碰撞啊 接触碰撞
		////如果算出的这条边包含原点 那么也肯定碰撞啊 


		//bool containsOrigin(simplex s, Vector d)
		//{
		//	//一般都是取刚刚添加的点为起点

		//	{
		//		//计算方向的时候还需要把离原点较远的点去掉
		//		//每次计算仅仅需要两个点

		//		//我们还要继续增加单纯形
		//		//d.getDirection(simplex);
		//		//选定第一个方向向量
		//		//Vector2 d;
		//		//Vector2 a = support(.., .., d);
		//		//Vector2 b = support(.., .., -d);
		//		//Vector2 AB = a - b;
		//		//Vector2 A0 = a - Vector2(0, 0);
		//		//d = AB X AO X AB; //这一句话我得自己算出来 它虽然是二维的 但是最终的结果也确实挥着到本来的平面
		//		// a x b x c = b(a*c) - a(b*c)
		//		//其实这里只是求了刚刚算出来的两个点连成的向量的一个垂直向量
		//		//难不成方向有讲究 这样算出来是最好的

		//		//!!精髓
		//		//这里的作用其实是通过两个点选取另外一个法向量
		//		//A B 的顺序是无所谓的
		//		//但是一定要满足这样的计算形式
		//		//BA X BO X BA
		//		//也可以
		//		//它会让下一次选择的向量指向原点那一方

		//	}


		//	//一定要获得刚刚添加的点 设为a
		//	//其余两点为b c 顺序无关
		//	//我们要判断的是 原点到底在 ab 外面 还是在 ac外面 还是在三角形里面
		//	//通过 acxabxab可以获得 ab 指向外面的法向量
		//	//通过 abxacxac可以获得 ac 指向外面的法向量
		//	//他们都有可能等于0啊 这个原点位于单纯行上到底意味这什么???
		//	//获得最近添加的点
		//	a = simplex.getLast();
		//	//compute AO
		//	ao = a.negate();
		//	if (simplex.points.size() == 3)
		//	{
		//		//then its the triangle case
		//		//get b and c
		//		b = simplex.getB();
		//		c = simplex.getC();
		//		//compute the edges
		//		ab = b - c;
		//		ac = c - a;
		//		//compute the normals
		//		abPrep = tripleProduct(ac, ab, ab);
		//		acPrep = tripleProduct(ab, ac, ac);
		//		//is the origin in R4
		//		if (abPrep.dot(ao) > 0)
		//		{
		//			//remove point c
		//			simplex.remove(c);
		//			//set the new direction to abPrep
		//			d.set(abPrep);
		//		}
		//		else
		//		{
		//			//is origin in R3
		//			if (acPrep.dot(ao) > 0)
		//			{
		//				//remove point b
		//				simplex.remove(b);
		//				//set the new direction to acPrep
		//				d.set(acPrep);
		//			}
		//			else
		//			{
		//				//otherwise we hnow its in R5 wo we can return 
		//				return true;
		//			}
		//		}
		//	}
		//	else
		//	{
		//		//then it is the line segment case
		//		b = simplex.getB();
		//		//compute AB
		//		ab = b - a;
		//		//get the perp to AB in the direction of the origin
		//		abPerp = tripleProduct(ab, ab, ab);

		//		//这里会出现一种情况
		//		//当通过两个点选择另外一个向量 另这个向量指向原点的方向时
		//		//这个向量可能是零向量
		//		//此时意味这 线段刚好通过原点
		//		/*
		//		The catch here is what happens when O lies on the line?
		//		If that happens the perp will be a zero vector and will cause the check on line 11 to fail.
		//		This can happen in two places: 1) inside the Minkowski Sum and 2) on the edge of the Minkowski Sum.
		//		The latter case indicates a touching contact rather than penetration so you will need to make a decision on whether this is considered a collision or not.
		//		In either case, you can use either the left or right hand normal of AB as the new direction.
		//		*/
		//		//有可能发生接触碰撞


		//		//set the direction to abPerp
		//		d.set(abPerp);
		//	}
		//}


		////---------------------------------计算距离和最近点-----------------------------------------

		////连个物体之间的最短距离就是它们的明可夫斯基差形状到原点的最近距离
		///*
		//如果两个凸体的明可夫斯基差形状中没有包括原点 则这两个物体不会相交
		//因此 我们不用采取迭代的方法不断计算包围原点的单纯形
		//而是想法得到一个离原点最近的单纯形
		//最近的单纯形总是在明可夫斯基差的边界上
		//*/

		///*
		//这一部分的代码看起来和GJK很想
		//不同的是 无论合适我们都只保留2个点 三维情况下是三个点
		//我们在单纯形中寻找离原点最近的点 而不是寻找原点所在Voronoi区域

		//这部分代码仅适用于凸多边形
		//而且 如果多边形相交 那么算法也会出错
		//不过一般没有问题 因为我们一般都会先检测碰撞

		//如果没有事先检测碰撞 我们就需要提供一个函数用来检测一个点是否在三角形内部

		//*/

		///*
		//另外 我们还可以确定两个多边形上最近的点
		//只不过我们需要额外存储一些信息
		//如果我们保存用来创建名科夫斯基差单纯形的两个多边形的点

		//声明四个变量
		//分别保存用来计算单纯形的 A B两点在原来的图形上的位置

		//我们就能用他们来确定两个多边形上最近的点
		//*/

		///* Convex Combination 凸组合
		//S是一个凸包
		//CH(S) = ∑λiPi = λ1P1 + ... + λnPn
		//where Pi∈S , λ∈R
		//and ∑λi = 1
		//where λi >= 0
		//*/

		///*
		//对于2D来说
		//CH(S) = λ1P1 + λ2P2 也就是说线段上的点可以用线段端点乘以特定系数来表示
		//如果Q是算法终止阶段单纯形上离原点最近的点
		//那么从Q到原点的向量一定垂直于Q所位于的线段
		//L = A - B
		//Q * L = 0

		//把Q的加权因子带入可以得到
		//(λ1*A + λ2 * B) * L = 0;
		//并且
		//λ1 + λ2 = 1
		//求解这个方程可以得到 λ1 λ2

		//求得加权因子之后 可以利用名科夫斯基差形状上的点来求得最近点
		//Acloset = λ1 * As1 + λ2 * Bs1  这里的两个点就是形成线段AB的图形S1上对应的点
		//Bcloset = λ1 * As2 + λ2 * Bs2
		//*/

		///*
		//仍有一些问题我们需要解决

		//1.如果明可夫斯基差形状上的A点和B点是一个点
		//那么L是零向量
		//这说明到原点的最近点不是在明可夫斯基形状的边上
		//而是其上的一个定点
		//这样A和B的support点是在同一个位置
		//我们可以返回A或B点
		//if(L.zero)
		//{
		//Acloset = A.s1;//A是一个support点 这个点由两个点计算而来 分别是A上的一个点 和 B上的一个点
		//Bcloset = A.s2;
		//}

		//2.第二个问题就是加权因子会出现负值的情况
		//如果为负 意味着另一个明可夫斯基差的support点为最近点
		//if(lambda1 < 0)
		//{
		//Acloest = B.s1;
		//Bcloest = B.s2;
		//}
		//else if(lambda2 < 0)
		//{
		//Acloset = A.s1;
		//Bcloset = A.s2;
		//}
		//*/

		////伪码
		//d = ;//choose a direction
		//	 //在单纯形中包含两个点
		//simplex.add(support(A, B, d));
		//simplex.add(support(A, B, -d));
		////开始迭代
		//while (true)
		//{
		//	//获得单纯形中两个点离原点最近的点
		//	p = closepointtoorigin(simplex.a, simplex.b);
		//	//检查是否为0
		//	if (p.zero())
		//	{
		//		//说明物体接触
		//		//物体接触意味这 距离为0  啊
		//		return false;
		//	}

		//	// p to origin 是新的向量
		//	//规范化是为了下面要用来投影计算
		//	d = p.negate().normallize();
		//	//在这个方向上获得一个新的名科夫斯基差点
		//	c = support(A, B, d);
		//	//新获得的点与之前放弃的那个点相比 我们更接近原点了吗?
		//	Real dc = c.dot(d);
		//	Real da = simpex.a.dot(d);
		//	//tolerance 是 精度 <= 0
		//	if (dc - da < tolerance)
		//	{
		//		//如果我们更接近了 
		//		//
		//		distance = dc;
		//		return true;
		//	}

		//	if (distance(a, ori) < distance(b, ori))
		//	{
		//		b = c;
		//	}
		//	else
		//	{
		//		a = c;
		//	}

		//}




		////----------------------------------穿透深度与穿透向量----------------------------------------

		///*
		//如果两个凸体的名科夫斯基差形状包括原点
		//则这两个物体是相交的
		//名科夫斯基差形状到原点的最短距离即为穿透深度
		//同样的
		//从最近点到原点形成的向量就是穿透向量
		//*/

		////EPA Expanding Polytope Algorigthm 扩张多面体算法
		///*
		//我们要在明可夫斯基差形状里面创建一个多面体 或多边形
		//迭代的扩张它 直到原点到单纯形最近距离的边就是明可夫斯基差形状的边
		//扩张的方法就是
		//根据多面体到原点的最短距离不多扩张多边形
		//通过迭代 我们产生的多面体包括明可夫斯基差形状到远点的最近点
		//这样就可以求得穿透深度和穿透向量
		//EPA的单纯形可以有任意多个点
		//*/
		//void EPA()
		//{
		//	/*
		//	启动点
		//	EPA算法在最开始的时候
		//	需要一个初始化的单纯形
		//	我们可以适用GJK算法终止时候的单纯形作为EPA算法的初始单纯形

		//	EPA算法需要完整的单纯形 在2d情况下是三角形
		//	在3d情况下是四面体
		//	*/
		//	simplex = ;//GJK算法终止的单纯形
		//	while (true)
		//	{
		//		//获得单纯形上距离原点最近的边
		//		Edge e = findcloseedge(s);
		//		//获得一个新的明可夫斯基差点 方向是这边的法向
		//		vector p = support(A, B, e.normal);

		//		//检查现在的这条边是不是就是离原点最近的边
		//		//如果差值为0 或者 非常小的话 我们都可以认为这就是那条边
		//		double d = p.dot(e.normal);
		//		if (d - e.distance < tolerance)
		//		{
		//			//tolerace 是浮点数精度 ex 1e-6
		//			//如果距离比精度还要小
		//			//那么我们可以认为我们不能把单纯形扩展的更远了
		//			//那么我们就获得了结果
		//			normal = e.normal;
		//			depth = d;
		//		}
		//		else
		//		{
		//			//我们没有找到
		//			//所以我们继续扩展单纯形
		//			//把这点添加到对应的下标位置 
		//			//simplex 大概要用链表这样的结构
		//			//为什么要添加到与下标相关的位置呢
		//			/*
		//			因为我们不应该破坏边的顺序
		//			就是说 我们在添加完这个点之后
		//			我们还是能按顺序来获取边
		//			*/
		//			simplex.insert(p, edge.index);
		//		}
		//	}
		//}

		////??这不就是遍历么...
		////这tm就是遍历找最小边
		///*

		//这个可以优化
		//用不着每次计算都重新计算所有的边
		//把已经计算过的边的值保存下来就行了
		//反正不缺那点内存
		//*/
		//void find_close_edge()
		//{
		//	edge closet;
		//	//把边的距离变成最大??
		//	closet.distance = Real_max_value;
		//	// s is the passed in simplex
		//	for (size_t i = 0; i < s.length; ++i)
		//	{
		//		//???
		//		size_t j = i + 1;
		//		j == s.length ? 0 : i + 1;//??
		//								  //get the current point and the next one
		//		vector a = s.get(i);
		//		vector b = s.get(j);
		//		//creat the edge vector
		//		vector e = b - a;
		//		//get vector OA
		//		vector oa = a - origin;
		//		//get the vector from edge towars the origin
		//		vector n = tripleProduct(e, oa, e);//计算远离原点那个方向的法向量
		//										   //normalize the vector
		//		n.mormalize();
		//		//calculate the distance from the origin to the dege
		//		// oooooooo  计算原点到这条边的距离
		//		double d = n.dot(a);
		//		//check the distance against the other distance
		//		if (d < closet.distance)
		//		{
		//			closet.distance = d;//原点距离最近边的距离
		//			closet.normal = n;//最近边的法向量
		//			closet.index = j;//第几条最近边
		//		}
		//	}
		//	//return the closet edge we found
		//	return closet;
		//}




		///*
		//需要注意的地方

		//算法中使用三重积得到一条边在背离原点方向的法向量
		//但是 在小的或者接触碰撞的情况下 使用三重积可能会使EPA算法产生问题
		//比如原点和最近边很近 三重积可能返回0向量

		//所以我们可以用边的垂直积来替代它
		//if A = (x,y)
		//A.preproduct() = (-y,x) or (y,-x)
		//依赖于坐标系

		//所以 不管原点离最近边多近 我们都能通过垂直积取得法向量
		//if(clockwise) (y,-x);
		//else (-y,x);

		//保持单纯形的旋转方向很重要
		//这样会使代码更有效

		//但是 什么使单纯形的旋转方向???
		//怎么保持???

		//...

		//没想到
		//但是我想到另外一个方向
		//仅仅需要执行一次向量乘法就可以得出结果
		//就是判断方向 然后在用 垂直积
		//*/


		///*
		//改进 ???

		//因为计算量很大 在穿透深度较小时 我们一般不适用EPA算法
		//所以EPA算法一般时GJK算法的补充
		//我们还可以用碰撞物体的一部分(核心形状)来进行检测和距离计算
		//一旦得到两个物体之间的距离 减去径向距离 就可以得到两个物体之间的距离
		//*

		////-------------------------------------------------------------------------------


} //namespace grid

#endif //_GRID_GEOMETRY_ALGORITHM_H_