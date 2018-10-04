//@header : gr_debug_draw.h
//@author : zhong
//@abstract : 用于简单的绘制各种几何图形

#ifndef _GRID_DEBUG_DRAW_H_
#define _GRID_DEBUG_DRAW_H_

#include<gr_geometry.h>
#include<graphics.h>
#include<string>

namespace grid
{
	namespace grid_debug
	{
		void DrawVector2(const Vector2& start, const Vector2& end)
		{
			if ((end - start).IsZero())
				return;
			line(
				static_cast<int>(start.x), static_cast<int>(start.y),
				static_cast<int>(end.x), static_cast<int>(end.y)
			);
			Vector2 vector = end - start;
			Real angle = 15.0;
			Real radian = AngleToRadian(angle);
			Real scale = 0.0678;
			Real length = vector.Length() * scale;
			Real height = length * GrTan(radian);
			Vector2 direction(vector.Invert().Normalize());
			Vector2 normal(VectorVerticalProduct(direction));
			Vector2 arrow[2]{};
			arrow[0] = end + direction * length + normal * height;
			arrow[1] = end + direction * length - normal * height;
			line(
				static_cast<int>(end.x), static_cast<int>(end.y),
				static_cast<int>(arrow[0].x), static_cast<int>(arrow[0].y)
			);
			line(
				static_cast<int>(end.x), static_cast<int>(end.y),
				static_cast<int>(arrow[1].x), static_cast<int>(arrow[1].y)
			);
		}

		void DrawPoint(const Point& point)
		{
			circle(
				static_cast<int>(point.x),
				static_cast<int>(point.y), 1);
		}

		void DrawSegment(const Segment& segment)
		{
			DrawVector2(segment.start, segment.end);
		}

		void DrawRay(const Ray& ray, Real distance)
		{
			Point end(ray.origin + ray.direction * distance);
			DrawVector2(ray.origin, end);
		}

		void DrawLine(const Line& line_,Real distance)
		{
			Vector2 direction(VectorVerticalProduct(line_.normal));
			Point start(line_.origin + direction * distance);
			Point end(line_.origin + direction * distance);
			line(static_cast<int>(start.x), static_cast<int>(start.y),
				static_cast<int>(end.x), static_cast<int>(end.y));
		}

		void DrawPolygon(const Polygon& poly)
		{
			POINT* points = new POINT[poly.n];
			if (points) {
				for (int i = 0; i < poly.n; ++i)
				{
					points[i].x = static_cast<LONG>(poly.vertex[i].x);
					points[i].y = static_cast<LONG>(poly.vertex[i].y);
				}
				polygon(points, poly.n);
				delete[] points;
			}
		}

		void DrawTriangle(const Triangle& triangle)
		{
			Polygon polytri(triangle);
			DrawPolygon(polytri);
		}

		void DrawRectangle(const Rectangle& rect)
		{
			rectangle(
				static_cast<int>(rect.x),
				static_cast<int>(rect.y),
				static_cast<int>(rect.x + rect.width),
				static_cast<int>(rect.y + rect.height)
			);
		}

		void DrawBox(const Box& box)
		{
			Polygon polybox(box);
			DrawPolygon(polybox);
		}

		void DrawArc(const Arc& arc)
		{
			Polygon polyarc(arc);
			DrawPolygon(polyarc);
		}

		void DrawSector(const Sector& sector)
		{
			Polygon polysec(sector);
			DrawPolygon(polysec);
		}

		void DrawCircle(const Circle& cir)
		{
			circle(
				static_cast<int>(cir.center.x),
				static_cast<int>(cir.center.y),
				static_cast<int>(cir.radius)
			);
		}

		void DrawString(Real x, Real y, const std::string& str)
		{
			outtextxy(
				static_cast<int>(x),
				static_cast<int>(y),
				str.c_str()
			);
		}

		void DrawString(const Vector2& position, const std::string& str)
		{
			DrawString(position.x, position.y, str);
		}

	} //namespace grid_debug

} //namespace grid

#endif //_GRID_DEBUG_DRAW_H_