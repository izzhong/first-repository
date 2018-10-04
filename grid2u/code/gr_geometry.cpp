//@source	:	gr_geometry.cpp
//@author	:	zhong
//@date		:	2018/9

#include"include/gr_geometry.h"
#include<iostream>

namespace grid
{
	PointArray::PointArray(const PointSet& pointset)
	{
		if (!pointset.Empty()){
			if (pointarray_ = Alloc_(pointset.Size())) {
				size_ = pointset.Size();
				PointNode* work = pointset.HeadPtr();
				for (int i = 0; i < size_; ++i,work = work->next)
				{
					pointarray_[i] = work->point;
				}
			}
		}
	}

	PointArray::PointArray(std::initializer_list<Point> il)
	{
		if (il.size() > 0)
		{
			if (pointarray_ = Alloc_(il.size())) {
				for (auto it = il.begin(); it != il.end(); ++it)
				{
					pointarray_[size_++] = *it;
				}
			}
		}		
	}

	PointArray::PointArray(const PointArray& pointarray)
	{
		if (!pointarray.Empty()) {
			if (pointarray_ = Alloc_(pointarray.Size())) {
				size_ = pointarray.size_;
				for (int i = 0;i < size_;++i)
				{
					pointarray_[i] = pointarray[i];
				}
			}
		}
	}

	//PointArray::PointArray(PointArray&& pointarray)
	//{
	//	if (!pointarray.Empty()) {			
	//		size_ = pointarray.size_;
	//		pointarray_ = pointarray_;
	//		pointarray.size_ = 0;
	//		pointarray.pointarray_ = nullptr;
	//	}
	//}

	PointArray& PointArray::operator=(const PointArray& pointarray)
	{
		Clear();
		if (!pointarray.Empty()) {
			if (pointarray_ = Alloc_(pointarray.Size())) {
				size_ = pointarray.size_;
				for (int i = 0; i < size_; ++i)
				{
					pointarray_[i] = pointarray[i];
				}
			}
		}
		return *this;
	}

	//PointArray& PointArray::operator=(PointArray&& pointarray)
	//{
	//	Clear();
	//	if (!pointarray.Empty()) {
	//		size_ = pointarray.size_;
	//		pointarray_ = pointarray_;
	//		pointarray.size_ = 0;
	//		pointarray.pointarray_ = nullptr;
	//	}
	//	return *this;
	//}

	PointSet::PointSet(const PointSet& pointset) :
		size_(0), head_(nullptr)
	{
		PointNode* work = pointset.head_;
		if (work)
		{
			for (int i = 0; i < pointset.Size(); ++i)
			{
				Push(work->point);
				work = work->next;
			}
		}
	}

	//PointSet::PointSet(PointSet&& pointset) :
	//	size_(0), head_(nullptr)
	//{
	//	head_ = pointset.head_;
	//	size_ = pointset.size_;
	//	pointset.head_ = nullptr;
	//	pointset.size_ = 0;
	//}

	PointSet::PointSet(int size, const Point* point_array) :
		size_(0), head_(nullptr)
	{
		Set(size, point_array);
	}

	PointSet& PointSet::Set(int size, const Point* point_array)
	{
		assert(size > 0);
		assert(point_array);
		Destruct_();
		for (int i = 0; i < size; ++i)
			Push(point_array[i]);
		return *this;
	}

	//PointSet& PointSet::operator=(const PointSet& pointset)
	//{
	//	Destruct_();
	//	PointNode* work = pointset.head_;
	//	if (work)
	//	{
	//		for (int i = 0; i < pointset.Size(); ++i)
	//		{
	//			Push(work->point);
	//			work = work->next;
	//		}
	//	}
	//	return *this;
	//}

	PointSet& PointSet::operator=(PointSet&& pointset)
	{
		Destruct_();
		head_ = pointset.head_;
		size_ = pointset.size_;
		pointset.head_ = nullptr;
		pointset.size_ = 0;
		return *this;
	}

	PointNode* PointSet::Alloc_(const Point& point)
	{
		return new PointNode{ point,nullptr,nullptr };
	}

	void PointSet::Destruct_()
	{
		while (!Empty())
		{
			EraseHead();
		}
	}

	void PointSet::Transform(const Matrix3& matrix)
	{
		if (Empty())
			return;
		PointNode* work = head_;
		do {
			matrix.Transform(&work->point);
			work = work->next;
		} while (work != head_);
	}

	void PointSet::ForEach(void(*todo)(Point*), bool positive)
	{
		if (Empty())
			return;
		PointNode* work = head_;
		do {
			todo(&work->point);
			if (positive) {
				work = work->next;
			}
			else {
				work = work->prev;
			}
		} while (work != head_);
	}

	void PointSet::ForEach(void(*todo)(const Point*), bool positive) const
	{
		if (Empty())
			return;
		PointNode* work = head_;
		do {
			todo(&work->point);
			if (positive) {
				work = work->next;
			}
			else {
				work = work->prev;
			}
		} while (work != head_);
	}

	void PointSet::ForEach(void(*todo)(Point&), bool positive)
	{
		if (Empty())
			return;
		PointNode* work = head_;
		do {
			todo(work->point);
			if (positive) {
				work = work->next;
			}
			else {
				work = work->prev;
			}
		} while (work != head_);
	}

	void PointSet::ForEach(void(*todo)(const Point&), bool positive) const
	{
		if (Empty())
			return;
		PointNode* work = head_;
		do {
			todo(work->point);
			if (positive) {
				work = work->next;
			}
			else {
				work = work->prev;
			}
		} while (work != head_);
	}

	PointNode* PointSet::Find(bool(*unary)(const Point&), bool positive)
	{
		if (Empty())
			return nullptr;
		PointNode* work = head_;
		do {
			if (unary(work->point))
				return work;
			if (positive) {
				work = work->next;
			}
			else {
				work = work->prev;
			}
		} while (work != head_);
		return nullptr;
	}

	PointNode* PointSet::Find(const Point& point, bool positive)
	{
		if (Empty())
			return nullptr;
		PointNode* work = head_;
		do {
			if (point == work->point)
				return work;
			if (positive) {
				work = work->next;
			}
			else {
				work = work->prev;
			}
		} while (work != head_);
		return nullptr;
	}

	PointNode* PointSet::InsertBefore(PointNode* at, const Point& point)
	{
		if (!at)
			return nullptr;
		PointNode* node = Alloc_(point);
		if (node)
		{
			node->next = at;
			node->prev = at->prev;
			at->prev->next = node;
			at->prev = node;
			if (at == head_)
				head_ = node;
			++size_;
		}
		return node;
	}

	PointNode* PointSet::InsertAfter(PointNode* at, const Point& point)
	{
		if (!at)
			return nullptr;
		PointNode* node = Alloc_(point);
		if (node)
		{
			node->next = at->next;
			node->prev = at;
			at->next->prev = node;
			at->next = node;
			++size_;
		}
		return node;
	}

	bool PointSet::Erase(const Point& point)
	{
		PointNode* work = Find(point);
		if (work) {
			Erase(work);
			return true;
		}
		else {
			return false;
		}
	}

	void PointSet::Erase(const PointNode* node)
	{
		if (Empty())
			return;
		if (node)
		{
			if (node == head_)
			{
				EraseHead();
				return;
			}
			node->prev->next = node->next;
			node->next->prev = node->prev;
			delete node;
			--size_;
		}
	}

	void PointSet::EraseTail()
	{
		if (Empty())
			return;
		if (Size() == 1)
		{
			delete head_;
			head_ = nullptr;
			--size_;
			return;
		}
		PointNode* tail = head_->prev;
		tail->prev->next = tail->next;
		tail->next->prev = tail->prev;
		delete tail;
		--size_;
	}

	void PointSet::EraseHead()
	{
		if (Empty())
			return;
		if (size_ == 1)
		{
			delete head_;
			head_ = nullptr;
			--size_;
			return;
		}
		PointNode* head = head_;
		head_ = head->next;
		head->prev->next = head_;
		head_->prev = head->prev;
		delete head;
		--size_;
	}

	PointNode* PointSet::Push(const Point& point)
	{
		PointNode* node = Alloc_(point);
		if (node)
		{
			if (Empty())
			{
				head_ = node;
				head_->next = head_;
				head_->prev = head_;
			}
			else
			{
				node->next = head_;
				node->prev = head_->prev;
				head_->prev->next = node;
				head_->prev = node;
			}
			++size_;
		}
		return node;
	}

	void PointSet::Push(PointNode* pointnode)
	{
		PointNode* node = pointnode;
		if (node)
		{
			if (Empty())
			{
				head_ = node;
				head_->next = head_;
				head_->prev = head_;
			}
			else
			{
				node->next = head_;
				node->prev = head_->prev;
				head_->prev->next = node;
				head_->prev = node;
			}
			++size_;
		}
	}

	PointSet& PointSet::Merge(const PointSet& pointset)
	{
		if (!pointset.Empty())
		{
			PointNode* work = pointset.head_;
			for (int i = 0; i < pointset.Size(); ++i)
			{
				Push(work->point);
				work = work->next;
			}
		}
		return *this;
	}

	PointSet& PointSet::Merge(PointSet&& pointset)
	{	
		if (!pointset.Empty()) 
		{
			PointNode* work = pointset.head_;
			int size = pointset.size_;
			for (int i = 0; i < size; ++i)
			{
				Push(work);
				work = work->next;
			}
			pointset.head_ = nullptr;
			pointset.size_ = 0;
		}	
		return *this;
	}

	Vector2 Segment::Normal(bool clock_wise) const
	{
		if (!IsStandard()) {
			return Vector2(0.0, 0.0);
		}
		Vector2 temp(end - start);
		temp.Normalize();
		if (clock_wise)
			return Vector2(temp.y, -temp.x);
		else
			return Vector2(-temp.y, temp.x);
	}

	Vector2 Triangle::SideNormal(int index) const
	{
		assert(index >= 0 && index < 3);
		Vector2 normal = VectorTripleProduct(
			vertex[(index + 1) % 3] - vertex[index],
			vertex[index] - vertex[(index + 2) % 3],
			vertex[(index + 1) % 3] - vertex[index]
		);
		normal.NormalizeSafe();
		return normal;
	}

	void Box::ToRectangle(Rectangle* rectangle, Matrix3* rotate_matrix) const
	{
		assert(rectangle);
		assert(rotate_matrix);
		Real radian = VectorRadian(direction_x, Vector2(1.0, 0.0));
		RotateMatrix(center, RadianToAngle(radian), rotate_matrix);
		rectangle->x = center.x - half_width;
		rectangle->y = center.y - half_height;
		rectangle->width = 2.0 * half_width;
		rectangle->height = 2.0 * half_height;
	}

	void Box::Transform(const Matrix3& matrix)
	{
		Point Dx(center + direction_x * half_width);
		Point Dy(center + direction_y * half_height);
		matrix.Transform(&center);
		matrix.Transform(&Dx);
		matrix.Transform(&Dy);
		half_width = (Dx - center).Length();
		half_height = (Dy - center).Length();
		direction_x = (Dx - center).NormalizeSafe();
		direction_y = (Dy - center).NormalizeSafe();
	}

	Box Rectangle::ToBox() const
	{
		return Box(
			Center(),
			Vector2(1, 0),
			Vector2(0, 1),
			width * 0.5,
			height * 0.5
		);
	}

	bool Polygon::IsStandard() const
	{
		Real wise =
			vertex[n - 1].x * vertex[0].y -
			vertex[n - 1].y * vertex[0].x;
		for (int i = 0; i < n - 1; ++i)
		{
			if (wise * (vertex[i].x * vertex[i + 1].y -
				vertex[i].y * vertex[i + 1].x) > 0.0) {
				continue;
			}
			else {
				return true;
			}
		}
		return true;
	}

	Point Polygon::Center() const
	{
		Point temp;
		for (int i = 0; i < n; ++i)
			temp += vertex[i];
		temp /= n;
		return temp;
	}

	Real Polygon::Area() const
	{
		Real area = vertex[n - 1].x * vertex[0].y - vertex[n - 1].y * vertex[0].x;
		for (int i = 0; i < n - 1; ++i)
		{
			area += vertex[i].x * vertex[i + 1].y - vertex[i].y * vertex[i + 1].x;
		}
		return std::fabs(0.5 * area);
	}

} //namespace grid