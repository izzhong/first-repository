//@header	:	gr_bounding.h
//@date		:	2018/9
//@author	:	zhong

//TIPS(zhong) : 不要滥用assert只需要在必要的时候提醒

#ifndef _GRID_BOUNDING_H_
#define _GRID_BOUNDING_H_

#include"gr_geometry_algorithm.h"

namespace grid
{
	class Bounding;
	class CenterBoundingCircle;
	class AxisAlignedBoundingBox;
	class OrientedBoundingBox;
	//class SpaceTimeBoundingBox;
	//class DisperseOrientedPolygon;
	using CBC = CenterBoundingCircle;
	using AABB = AxisAlignedBoundingBox;
	using OBB = OrientedBoundingBox;
	//using kDOP = DisperseOrientedPolygon;
	//using STBB = SpaceTimeBoundingBox;

	//@abstract : Bounding
	class Bounding
	{	
	public:
		virtual void CreateBounding(
			int size, const Point* pointarray
		) = 0;

		virtual void CreateBounding(
			const PointSet& pointset
		) = 0;

		virtual void CreateBounding(
			const Circle& circle
		) = 0;

		virtual void CreateBounding(
			const Arc& arc
		) = 0;

		virtual void CreateBounding(
			const Sector& sector
		) = 0;

		virtual void CreateBounding(
			const Triangle& triangle
		) = 0;

		virtual void CreateBounding(
			const Rectangle& rectangle
		) = 0;

		virtual void CreateBounding(
			const Box& box
		) = 0;

		virtual void CreateBounding(
			const Polygon& polygon
		) = 0;

		virtual void Clear() = 0;

		virtual void Transform(const Matrix3& matrix) = 0;
	};

	//@abstract : cbc
	class CenterBoundingCircle :
		public Bounding
	{
	public:

		CenterBoundingCircle() {  }

		CenterBoundingCircle(const Circle& circle) :
			circle_(circle) {	}

		virtual void Clear() override
		{
			circle_.Clear();
		}

		Circle GetGeometry() const
		{
			return Circle(circle_);
		}

		virtual void Transform(const Matrix3& matrix) override
		{
			circle_.Transform(matrix);
		}

		bool Seperate(const CenterBoundingCircle& cbc) const
		{
			return CircleSeparationDetection(
				circle_,cbc.GetGeometry()
			);
		}

		friend CenterBoundingCircle CombineBounding(
			const CenterBoundingCircle& cbc_1,
			const CenterBoundingCircle& cbc_2
		);

		//bool Seperate(const AxisAlignedBoundingBox& aabb) const;

		//bool Seperate(const OrientedBoundingBox& obb) const;

	public:

		virtual void CreateBounding(
			int size, const Point* point_array
		) override;
		
		virtual void CreateBounding(
			const PointSet& pointset
		) override;

		virtual void CreateBounding(
			const Circle& circle
		) override;

		virtual void CreateBounding(
			const Triangle& triangle
		) override;

		virtual void CreateBounding(
			const Rectangle& rectangle
		) override;

		virtual void CreateBounding(
			const Box& box
		) override;

		virtual void CreateBounding(
			const Polygon& polygon
		) override;

		virtual void CreateBounding(
			const Arc& arc
		) override;

		virtual void CreateBounding(
			const Sector& sector
		) override;

	private:

		Circle circle_;

	private:

		void CreateCenterBoundingCircle(
			int size, const Point* point_array,
			Circle* cb
		);
	};

	class AxisAlignedBoundingBox :
		public Bounding
	{
	public:

		AxisAlignedBoundingBox() {	}
	
		AxisAlignedBoundingBox(const Rectangle& rectangle) :
			rect_(rectangle) {	}

		virtual void Clear() override
		{
			rect_.Clear();
		}
	
		Rectangle GetGeometry() const
		{
			return rect_;
		}

		bool Seperate(const AxisAlignedBoundingBox& bounding)
		{
			//直接调用几何学算法
			return RectangleSeparationDetection(
				rect_, bounding.GetGeometry()
			);
		}

		virtual void Transform(const Matrix3& matrix) override
		{
			rect_.Transform(matrix);
		}

		//bool Seperate(const CBC& cbc);

		//bool Seperate(const OBB& obb);

		friend AxisAlignedBoundingBox CombineBounding(
			const AxisAlignedBoundingBox& aabb_1,
			const AxisAlignedBoundingBox& aabb_2
		);

	public:

		virtual void CreateBounding(
			int size, const Point* point_array
		) override;

		virtual void CreateBounding(
			const PointSet& point_set
		) override;

		virtual void CreateBounding(
			const Circle& circle
		) override;

		virtual void CreateBounding(
			const Triangle& triangle
		) override;

		virtual void CreateBounding(
			const Rectangle& rectangle
		) override;

		virtual void CreateBounding(
			const Box& box
		) override;

		virtual void CreateBounding(
			const Polygon& polygon
		) override;

		virtual void CreateBounding(
			const Arc& arc
		) override;

		virtual void CreateBounding(
			const Sector& sector
		) override;

	private:

		Rectangle rect_;
	};

	class OrientedBoundingBox :
		public Bounding
	{
	public:

		OrientedBoundingBox() { }

		OrientedBoundingBox(const Box& box) :
			box_(box) {	}

		Box GetGeometry() const
		{
			return box_;
		}

		virtual void Clear() override
		{
			box_.Clear();
		}

		virtual void Transform(const Matrix3& matrix) override
		{
			box_.Transform(matrix);
		}

		friend OrientedBoundingBox CombineBounding(
			const OrientedBoundingBox& obb_1,
			const OrientedBoundingBox& obb_2
		);

		virtual bool Seperate(const OrientedBoundingBox& obb)
		{
			return false;
		}

	public:

		virtual void CreateBounding(
			int size, const Point* pointarray
		) override;

		virtual void CreateBounding(
			const PointSet& pointset
		) override;

		virtual void CreateBounding(
			const Circle& circle
		) override;

		virtual void CreateBounding(
			const Arc& arc
		) override;

		virtual void CreateBounding(
			const Sector& sector
		) override;

		virtual void CreateBounding(
			const Triangle& triangle
		) override;

		virtual void CreateBounding(
			const Rectangle& rectangle
		) override;

		virtual void CreateBounding(
			const Box& box
		) override;

		virtual void CreateBounding(
			const Polygon& polygon
		) override;

	private:

		Box box_;
	};

	//创建 k-DOP包围体的方法
	//该方法可以推广到所有凸多边形包围体的创建
	//将顶点投影到凸多边形的每条法线方向ni上，
	//然后把投影的极值保存在dimin,dimax中
	//这两个值定义了这个方向上最紧密的平板层 夹层
	//这些值共同定义了一个最小的凸多边形

} //namespace grid

#endif //_GRID_BOUNDING_H_