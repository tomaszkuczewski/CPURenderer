#pragma once

inline float clampf(float min, float max, float value)
{
	if (value < min) return min;
	else if (value > max) return max;
	else return value;
}

struct Mat44;
struct Vec4;

struct Vec3
{
	union
	{
		struct { float x, y, z; };
		struct { float _m[3]; };
	};

	inline Vec3 operator*(float value) const;
	inline Vec3 operator/(float value) const;
	inline Vec3 operator+(const Vec3& vector3) const;
	inline Vec3 operator-(const Vec3& vector3) const;

	inline Vec3& operator*=(float value);
	inline Vec3& operator/=(float value);
	inline Vec3& operator+=(const Vec3& vector3);
	inline Vec3& operator-=(const Vec3& vector3);

	inline void ToVec4(Vec4* vector, float w);
	inline float Dot(const Vec3& vector3) const;
	inline Vec3 Cross(const Vec3& vector3) const;
	inline void Normalize();
	inline float Length();

	inline static bool IntersectWithPlane(Vec3& out, float& tOut, const Vec3& planeNormal, const Vec3& planePoint, const Vec3& rayPoint, const Vec3& rayDir);
};

struct Vec4
{
	union
	{
		struct { float x, y, z, w; };
		struct { float _m[4]; };
	};

	inline Vec4 operator*(float value) const;
	inline Vec4 operator/(float value) const;
	inline Vec4 operator+(const Vec4& vector3) const;
	inline Vec4 operator-(const Vec4& vector3) const;

	inline Vec4& operator*=(float value);
	inline Vec4& operator/=(float value);
	inline Vec4& operator+=(const Vec4& vector3);
	inline Vec4& operator-=(const Vec4& vector3);

	//For matrix operations
	inline Vec4 operator*(const Mat44& mat) const;
	inline Vec4& operator*=(const Mat44& mat);

	//NOTE: Added only to make easier to write
	//multiplication of 4x4 matrices or vectors * matrices
	inline float Dot(const Vec4& vector4) const;

	inline void ToVec3(Vec3& output) const;
};

struct Mat44
{
	float _m[4][4];

	//For easier calculations
	inline void ToColumn(Vec4& output, int column) const;
	inline void ToRow(Vec4& output, int row) const;

	//Views
	inline static Mat44 ViewLH(const Vec3& eye, const Vec3& target, const Vec3& up);

	//Projections
	inline static Mat44 ProjectionLH(float width, float height, float fov, float zNear, float zFar);

	//Rotations
	inline static Mat44 RotationY(float value);
	inline static Mat44 RotationX(float value);
	inline static Mat44 RotationZ(float value);

	//For easier calculations
	Mat44& operator*=(const Mat44& mat);
	Mat44 operator*(const Mat44& mat) const;
};

#include <chrono>
#include <vector>
#include <functional>

class Mesh
{
private:
	std::vector<int> _mIndices;
	std::vector<float> _mVertices;
public:
	Mesh() = default;
	virtual ~Mesh() = default;

	inline std::vector<int>& GetIndices() { return _mIndices; }
	inline std::vector<float>& GetVertices() { return _mVertices; }
};

class CPURenderer
{
private:
	bool _mIsWorking = false;

	uint32_t _mWidth;
	uint32_t _mHeight;

	Vec3 _mPlaneNormals[6];
	Vec3 _mPlanePoints[6];

	const Mat44* _mLocalWorld;
	const Mat44* _mLocalView;
	const Mat44* _mLocalProjection;

	bool _mIsLightSet;
	Vec3 _mLight;

	//For performance calculations
	std::chrono::high_resolution_clock::time_point _mTimeNow;
	float _mfDeltaTime;
	float _mfFPS;
	
	std::vector<uint32_t> _mvBackBuffer;

	std::function<void(float)> _mFOnUpdate;
	std::function<void(uint32_t, uint32_t, uint32_t)> _mFOnDrawPixel;

	inline void SamplePerformance();
	inline void DrawPixelBuffered(uint32_t x, uint32_t y, uint32_t color);
	inline void DrawLineBuffered(int x1, int y1, int x2, int y2, uint32_t color);

	inline void DrawPixels();

protected:
	virtual int ClipInProjection(const Vec3* v_in, int inputTriangles, Vec3* v_out) const;
	virtual int ClipAgainstPlane(const Vec3* v, const Vec3* planeNormal, const Vec3* planePoint, Vec3* newTriangle0) const;
	virtual void DrawTriangleFrame(const Vec3* v0, const Vec3* v1, const Vec3* v2, uint32_t color);
	virtual void DrawTriangle(const Vec3* v0, const Vec3* v1, const Vec3* v2, uint32_t color);
	virtual void FillFlatBottomTriangle(const Vec3* v0, const Vec3* v1, const Vec3* v2, uint32_t color);
	virtual void FillFlatTopTriangle(const Vec3* v0, const Vec3* v1, const Vec3* v2, uint32_t color);

	inline void SortVerticesByY (Vec3* vertices, int count) const;
public:
	CPURenderer();
	CPURenderer(uint32_t width, uint32_t height);
	virtual ~CPURenderer() = default;

	inline int GetWidth() const noexcept { return _mWidth; }
	inline int GetHeight() const noexcept  { return _mHeight; }
	inline float GetWidthf() const noexcept { return static_cast<float>(_mWidth); }
	inline float GetHeightf() const noexcept  { return static_cast<float>(_mHeight); }

	inline void Clear();
	inline void Draw(Mesh* mesh);

	inline void Start			(const std::function<bool(int, int)>& startFunction);
	inline void Stop() { _mIsWorking = false; }

	//Seters for matrices
	inline void SetWorldMatrix(const Mat44* world) { _mLocalWorld = world; }
	inline void SetViewMatrix(const Mat44* world) { _mLocalView = world; }
	inline void SetProjectionMatrix(const Mat44* world) { _mLocalProjection = world; }

	//Setters for lights
	inline void SetLight(Vec3 lightDir) { _mLight = lightDir; _mIsLightSet = true; }

	//Setters for pixels
	inline void SetDrawPixel	(const std::function<void(uint32_t, uint32_t, uint32_t)>& drawPixelFunction) { _mFOnDrawPixel = drawPixelFunction; }
	inline void SetUpdate		(const std::function<void(float)>& updateFunction) { _mFOnUpdate = updateFunction; }
};

inline CPURenderer::CPURenderer() : CPURenderer{ 640, 480 }
{
}

inline CPURenderer::CPURenderer(uint32_t width, uint32_t height) :
	_mWidth{ width }, _mHeight{ height },
	_mLocalWorld{ nullptr }, _mLocalView{ nullptr }, _mLocalProjection{ nullptr },
	_mfDeltaTime{ 0.0f }, _mfFPS{0.0f}
{
	_mPlaneNormals[0] = Vec3{ 0, 0, -1 }; _mPlanePoints[0] = Vec3{ 0, 0, 1 }; 	//Far Z
	_mPlaneNormals[1] = Vec3{ 1, 0, 0 }; _mPlanePoints[1] = Vec3{ 0, 0, 0 }; //Left
	_mPlaneNormals[2] = Vec3{ -1, 0, 0 }; _mPlanePoints[2] = Vec3{ GetWidthf(), 0, 0 }; //Right

	_mvBackBuffer.resize(width * height);
	_mvBackBuffer.shrink_to_fit();
}

inline void CPURenderer::Start(const std::function<bool(int, int)>& startFunction)
{
	_mIsWorking = true;

	if (!_mFOnUpdate || !_mFOnDrawPixel || !startFunction )
		return;

	//Call the start lamda
	if (!startFunction(GetWidth(), GetHeight()))
		return;

	while (_mIsWorking)
	{
		SamplePerformance(); //Getting the FPS
		_mFOnUpdate(_mfDeltaTime);

		DrawPixels();
	}
}

inline void CPURenderer::Clear()
{
	std::fill(_mvBackBuffer.begin(), _mvBackBuffer.end(), 0);
}

inline void CPURenderer::Draw(Mesh* mesh)
{
	static size_t strideSize = 7;

	if (_mLocalView == nullptr || _mLocalProjection == nullptr)
		return;

	Vec4 vertex[6]; //Preallocated vertices
	Vec3 vertexClipped[24]; //Vertices clipped by view planes
	Vec3 vertexColor{};

	//Nera plane clip plane to remove artifacts
	int nearPlanceClippedCount = 0;
	auto nearPlanePoint = Vec3{ 0, 0, 0.1f };
	auto nearPlaneNormal = Vec3{ 0, 0, 1 };

	//Foreach index
	for (size_t i = 0; i < mesh->GetIndices().size(); i+=3)
	{
		//Input merger v0
		vertex[0].x = mesh->GetVertices()[ mesh->GetIndices()[i] * strideSize + 0 ];
		vertex[0].y = mesh->GetVertices()[ mesh->GetIndices()[i] * strideSize + 1 ];
		vertex[0].z = mesh->GetVertices()[ mesh->GetIndices()[i] * strideSize + 2 ];
		vertex[0].w = 1.0f;

		//Input merger v1
		vertex[1].x = mesh->GetVertices()[mesh->GetIndices()[i+1] * strideSize + 0];
		vertex[1].y = mesh->GetVertices()[mesh->GetIndices()[i+1] * strideSize + 1];
		vertex[1].z = mesh->GetVertices()[mesh->GetIndices()[i+1] * strideSize + 2];
		vertex[1].w = 1.0f;

		//Input merger v2
		vertex[2].x = mesh->GetVertices()[mesh->GetIndices()[i+2] * strideSize + 0];
		vertex[2].y = mesh->GetVertices()[mesh->GetIndices()[i+2] * strideSize + 1];
		vertex[2].z = mesh->GetVertices()[mesh->GetIndices()[i+2] * strideSize + 2];
		vertex[2].w = 1.0f;

		//Vertex color
		vertexColor.x = mesh->GetVertices()[mesh->GetIndices()[i] * strideSize + 3];
		vertexColor.y = mesh->GetVertices()[mesh->GetIndices()[i] * strideSize + 4];
		vertexColor.z = mesh->GetVertices()[mesh->GetIndices()[i] * strideSize + 5];

		//To World transformation
		if (_mLocalWorld != nullptr)
		{
			vertex[0] *= *_mLocalWorld; 
			vertex[1] *= *_mLocalWorld;
			vertex[2] *= *_mLocalWorld;
		}

		//ONLY FOR FLAT COLOR
		//Handle lighting and sience is flat color lighting
		if (_mIsLightSet)
		{
			//Aka Pixel/Fragment shader
			Vec3 v01, v02;
			(vertex[1] - vertex[0]).ToVec3(v01);
			(vertex[2] - vertex[0]).ToVec3(v02);
			auto normal = v01.Cross(v02);
			normal.Normalize();
			vertex[0].ToVec3(v01);

			if (normal.Dot(v01) <= 0)
				continue;

			_mLight.Normalize();
			auto dot = clampf(0.0f, 1.0f, normal.Dot(_mLight));
			vertexColor = (vertexColor * dot);
		}

		//To View transformation
		vertex[0] *= (*_mLocalView);
		vertex[1] *= (*_mLocalView);
		vertex[2] *= (*_mLocalView);

		//Shrink vertices to vector3 for easier calculation
		vertex[0].ToVec3(vertexClipped[0]);
		vertex[1].ToVec3(vertexClipped[1]);
		vertex[2].ToVec3(vertexClipped[2]);

		//Get the clipped triangles by near plane
		if (!(nearPlanceClippedCount = ClipAgainstPlane(vertexClipped, &nearPlaneNormal, &nearPlanePoint, vertexClipped)))
			continue;

		//Go back to Vec4
		for (int i = 0; i < nearPlanceClippedCount; i++)
		{
			int idx0 = i * 3;
			int idx1 = i * 3 + 1;
			int idx2 = i * 3 + 2;
			vertexClipped[idx0].ToVec4(&vertex[idx0], 1.0f);
			vertexClipped[idx1].ToVec4(&vertex[idx1], 1.0f);
			vertexClipped[idx2].ToVec4(&vertex[idx2], 1.0f);

			//To Projection transformation
			vertex[idx0] *= *_mLocalProjection; vertex[idx1] *= *_mLocalProjection; vertex[idx2] *= *_mLocalProjection;

			//To NDC
			vertex[idx0] /= vertex[idx0].w;		vertex[idx1] /= vertex[idx1].w;					vertex[idx2] /= vertex[idx2].w; //Divide by w

			vertex[idx0].x += 1.0f;				vertex[idx0].x *= 0.5f;		vertex[idx0].y += 1.0f;		vertex[idx0].y *= 0.5f;
			vertex[idx0].x *= GetWidthf();		vertex[idx0].y *= GetHeightf();

			vertex[idx1].x += 1.0f;			vertex[idx1].x *= 0.5f;	vertex[idx1].y += 1.0f;	vertex[idx1].y *= 0.5f;
			vertex[idx1].x *= GetWidthf();	vertex[idx1].y *= GetHeightf();

			vertex[idx2].x += 1.0f;			vertex[idx2].x *= 0.5f;	vertex[idx2].y += 1.0f;	vertex[idx2].y *= 0.5f;
			vertex[idx2].x *= GetWidthf();	vertex[idx2].y *= GetHeightf();

			//Switch back to Vec3
			vertex[idx0].ToVec3(vertexClipped[idx0]); vertex[i * 3 + 1].ToVec3(vertexClipped[i * 3 + 1]); vertex[i * 3 + 2].ToVec3(vertexClipped[i * 3 + 2]);


			auto finalColor = (((uint32_t)(255.0f * vertexColor.x)) << 24) |
				(((uint32_t)(255.0f * vertexColor.y)) << 16) |
				(((uint32_t)(255.0f * vertexColor.z)) << 8) | 0xFF;

			SortVerticesByY(&vertexClipped[idx0], 3);
			DrawTriangle(&vertexClipped[idx0], &vertexClipped[idx1], &vertexClipped[idx2], finalColor);
		}
	}
}

inline void CPURenderer::SamplePerformance()
{
	//Calculate FPS
	_mfDeltaTime = (float)std::chrono::duration_cast<std::chrono::microseconds>
		(std::chrono::high_resolution_clock::now() - _mTimeNow).count() / 1000000.0f;
	//And then calulate fps 
	_mfFPS = 1.0f / _mfDeltaTime;
	//For performance
	_mTimeNow = std::chrono::high_resolution_clock::now(); //Performance measurements
}

inline void CPURenderer::DrawPixelBuffered(uint32_t x, uint32_t y, uint32_t color)
{
	_mvBackBuffer[ (y * _mWidth) + x] = color;
}

inline void CPURenderer::DrawLineBuffered(int x1, int y1, int x2, int y2, uint32_t color)
{
	x1 = (int)clampf(0, 1.0f* GetWidth()-1, x1 * 1.0f);
	x2 = (int)clampf(0, 1.0f* GetWidth()-1, x2 * 1.0f);
	y1 = (int)clampf(0, 1.0f * GetHeight()-1, y1 * 1.0f);
	y2 = (int)clampf(0, 1.0f * GetHeight()-1, y2 * 1.0f);

	if (x1 > x2)
	{
		std::swap(x2, x1);
		std::swap(y2, y1);
	}

	float dx = (1.0f * x2 - x1);
	float coeffA = (dx == 0.0f ? 0.0f : (1.0f * y2 - y1) / dx);
	for (int i = x1; i <= x2; i++)
		DrawPixelBuffered(i, y1 + (int)(coeffA * i), color);
}

inline void CPURenderer::DrawPixels()
{
	for (uint32_t x = 0; x < _mWidth; x++)
		for (uint32_t y = 0; y < _mHeight; y++)
			_mFOnDrawPixel(x, _mHeight - y, _mvBackBuffer[(y * _mWidth) + x]);
}

inline int CPURenderer::ClipInProjection(const Vec3* v_in, int inputTriangles, Vec3* v_out) const
{
	Vec3 calculationBuffer[24];

	int outTriangleCount = inputTriangles;
	int newTriangleCount;

	for (int i = 0; i < inputTriangles * 3; i++)
		v_out[i] = v_in[i];

	for (int i = 0; i < 3; i++)
	{
		newTriangleCount = 0;
		for (int t = 0; t < outTriangleCount; t++)
		{
			int ret = ClipAgainstPlane(v_out + ((long)t * 3), &_mPlaneNormals[i], &_mPlanePoints[i], (calculationBuffer + (3 * newTriangleCount)));
			if (ret == -1) return 0;
			else newTriangleCount += ret;
		}
		//Copy new triangles to out buffer
		for (int n = 0; n < newTriangleCount * 3; n++)
			v_out[n] = calculationBuffer[n];
		outTriangleCount = newTriangleCount;
	}
	return outTriangleCount;
}

inline int CPURenderer::ClipAgainstPlane(const Vec3* v, const Vec3* planeNormal, const Vec3* planePoint, Vec3* newTriangle0) const
{
	Vec3 v_copy[3];
	v_copy[0] = v[0];
	v_copy[1] = v[1];
	v_copy[2] = v[2];

	float t0, t1;
	Vec3 intersection0, intersection1;
	int insideCount = 0, outsideCount = 0;
	int insideVertices[3] = { -1, -1, -1 };
	int outsideVertices[3] = { -1, -1, -1 };

	//Check which vertices are outside/inside
	if (planeNormal->Dot(v[0] - *planePoint) <= 0) { outsideVertices[outsideCount++] = 0; }
	else insideVertices[insideCount++] = 0;
	if (planeNormal->Dot(v[1] - *planePoint) <= 0) { outsideVertices[outsideCount++] = 1; }
	else insideVertices[insideCount++] = 1;
	if (planeNormal->Dot(v[2] - *planePoint) <= 0) { outsideVertices[outsideCount++] = 2; }
	else insideVertices[insideCount++] = 2;

	//Render whole triangle
	if (outsideCount == 0)
	{
		newTriangle0[0] = v_copy[0]; newTriangle0[1] = v_copy[1]; newTriangle0[2] = v_copy[2];
		return 1;
	}
	//Skip whole triangle
	if (outsideCount == 3)
		return -1;
	//Make two new triangles
	if (outsideCount == 1)
	{
		auto toV0 = v_copy[outsideVertices[0]] - v_copy[insideVertices[0]]; toV0.Normalize();
		auto toV1 = v_copy[outsideVertices[0]] - v_copy[insideVertices[1]]; toV1.Normalize();
		Vec3::IntersectWithPlane(intersection0, t0, *planeNormal, *planePoint, v_copy[insideVertices[0]], toV0);
		Vec3::IntersectWithPlane(intersection1, t1, *planeNormal, *planePoint, v_copy[insideVertices[1]], toV1);
		newTriangle0[1] = v_copy[insideVertices[0]]; newTriangle0[0] = intersection0; newTriangle0[2] = v_copy[insideVertices[1]];
		newTriangle0[3] = intersection1; newTriangle0[4] = intersection0; newTriangle0[5] = v_copy[insideVertices[1]];
		return 2;
	}
	//Make one new triangle
	if (outsideCount == 2)
	{
		auto toV0 = v_copy[outsideVertices[0]] - v_copy[insideVertices[0]]; toV0.Normalize();
		auto toV1 = v_copy[outsideVertices[1]] - v_copy[insideVertices[0]]; toV1.Normalize();
		Vec3::IntersectWithPlane(intersection0, t0, *planeNormal, *planePoint, v_copy[insideVertices[0]], toV0);
		Vec3::IntersectWithPlane(intersection1, t1, *planeNormal, *planePoint, v_copy[insideVertices[0]], toV1);
		newTriangle0[0] = v_copy[insideVertices[0]]; newTriangle0[1] = intersection0; newTriangle0[2] = intersection1;
		return 1;
	}
	return -1;
}

inline void CPURenderer::DrawTriangleFrame(const Vec3* v0, const Vec3* v1, const Vec3* v2, uint32_t color)
{
	DrawLineBuffered((uint32_t)v0->x, (uint32_t)v0->y, (uint32_t)v1->x, (uint32_t)v1->y, color);
	DrawLineBuffered((uint32_t)v1->x, (uint32_t)v1->y, (uint32_t)v2->x, (uint32_t)v2->y, color);
	DrawLineBuffered((uint32_t)v2->x, (uint32_t)v2->y, (uint32_t)v0->x, (uint32_t)v0->y, color);
}

inline void CPURenderer::DrawTriangle(const Vec3* v0, const Vec3* v1, const Vec3* v2, uint32_t color)
{
	if (v1->y == v2->y)
		FillFlatBottomTriangle(v0, v1, v2, color);
	else if (v0->y == v1->y)
		FillFlatTopTriangle(v0, v1, v2, color);
	else
	{
		Vec3 v4{};
		float dx02 = (v2->x - v0->x);

		if (dx02 == 0.0)
		{
			v4.x = v0->x;
			v4.y = v1->y;
		}
		else
		{
			float a02 = (v2->y - v0->y) / dx02; 
			float b02 = v0->y - (a02)* v0->x;
			v4.x = (v1->y - b02) / a02;
			v4.y = (v1->y);
		}

		FillFlatBottomTriangle(v0, v1, &v4, color);
		FillFlatTopTriangle(&v4, v1, v2, color);
	}
}

inline void CPURenderer::FillFlatBottomTriangle(const Vec3* v0, const Vec3* v1, const Vec3* v2, uint32_t color)
{
	float invA01 = (v1->x - v0->x) / (v1->y - v0->y); //Inverse linear fnc coeff 'a' from 0->1
	float invA02 = (v2->x - v0->x) / (v2->y - v0->y); //Inverse linear fnc coeff 'a' from 0->2
	float currX0 = v0->x;
	float currX1 = v0->x;
	for (uint32_t y = (uint32_t) v0->y; y >= (uint32_t)v2->y; y--) //Scanline from v0->y to v2.y by decreasing y
	{
		DrawLineBuffered((uint32_t)currX0, y, (uint32_t)currX1, y, color);
		currX0 -= invA02;//Apply coeffs with change of y
		currX1 -= invA01;
	}
}

inline void CPURenderer::FillFlatTopTriangle(const Vec3* v0, const Vec3* v1, const Vec3* v2, uint32_t color)
{
	float invA01 = (v0->x - v2->x) / (v0->y - v2->y);
	float invA02 = (v1->x - v2->x) / (v1->y - v2->y);
	float currX0 = v2->x;
	float currX1 = v2->x;
	for (uint32_t y = (uint32_t)v2->y; y <= (uint32_t) v0->y; y++)
	{
		DrawLineBuffered((uint32_t)currX0, (uint32_t)y, (uint32_t)currX1, (uint32_t)y, color);
		currX0 += invA02;//Apply coeffs with change of y
		currX1 += invA01;
	}
}

inline void CPURenderer::SortVerticesByY(Vec3* vertices, int count) const
{
	std::qsort((void*)vertices, count, sizeof Vec3, [](const void* a, const void* b)
	{
		if (((Vec3*)a)->y > ((Vec3*)b)->y)
			return -1;
		else
			return 1;
	});
}

//Definitions
//#define _USE_MATH_DEFINES
#include <math.h>

inline Vec3 Vec3::operator*(float value) const
{
	return Vec3{ _m[0] * value, _m[1] * value , _m[2] * value };
}

inline Vec3 Vec3::operator/(float value) const
{
	return Vec3{ _m[0] / value, _m[1] / value , _m[2] / value };
}

inline Vec3 Vec3::operator+(const Vec3& vector3) const
{
	return Vec3{ _m[0] + vector3._m[0], _m[1] + vector3._m[1], _m[2] + vector3._m[2] };
}

inline Vec3 Vec3::operator-(const Vec3& vector3) const
{
	return Vec3{ _m[0] - vector3._m[0], _m[1] - vector3._m[1], _m[2] - vector3._m[2] };
}

inline Vec3& Vec3::operator*=(float value)
{
	_m[0] *= value; _m[1] *= value; _m[2] *= value;
	return *this;
}

inline Vec3& Vec3::operator/=(float value)
{
	_m[0] /= value; _m[1] /= value; _m[2] /= value;
	return *this;
}

inline Vec3& Vec3::operator+=(const Vec3& vector3)
{
	_m[0] += vector3._m[0]; _m[1] += vector3._m[1]; _m[2] += vector3._m[2];
	return *this;
}

inline Vec3& Vec3::operator-=(const Vec3& vector3)
{
	_m[0] -= vector3._m[0]; _m[1] -= vector3._m[1]; _m[2] -= vector3._m[2];
	return *this;
}

inline void Vec3::ToVec4(Vec4* vector, float w)
{
	*((Vec3*)vector) = *this;
	vector->w = w;
}

inline float Vec3::Dot(const Vec3& vector3) const
{
	return (_m[0] * vector3._m[0]) + (_m[1] * vector3._m[1]) + (_m[2] * vector3._m[2]);
}

inline Vec3 Vec3::Cross(const Vec3& vector3) const
{
	return Vec3{
			_m[1] * vector3._m[2] - _m[2] * vector3._m[1],
			_m[2] * vector3._m[0] - _m[0] * vector3._m[2],
			_m[0] * vector3._m[1] - _m[1] * vector3._m[0] };
}

inline void Vec3::Normalize()
{
	float d = sqrtf((_m[0] * _m[0]) + (_m[1] * _m[1]) + (_m[2] * _m[2]));
	this->operator/=(d);
}

inline float Vec3::Length()
{
	return sqrtf((_m[0] * _m[0]) + (_m[1] * _m[1]) + (_m[2] * _m[2]));
}

inline bool Vec3::IntersectWithPlane(Vec3& out, float& tOut, const Vec3& planeNormal, const Vec3& planePoint, const Vec3& rayPoint, const Vec3& rayDir)
{
	float denominator = rayDir.Dot(planeNormal);
	if (denominator == 0.0f)
		return false;

	float t = -(rayPoint.Dot(planeNormal) + (-planeNormal.Dot(planePoint))) / denominator;
	if (t == 0.0) return false;
	out = rayPoint + (rayDir * t);
	tOut = t;
	return true;
}

inline Vec4 Vec4::operator*(float value) const
{
	return Vec4{ _m[0] * value, _m[1] * value , _m[2] * value, _m[3] * value };
}

inline Vec4 Vec4::operator/(float value) const
{
	return Vec4{ _m[0] / value, _m[1] / value , _m[2] / value, _m[3] / value };
}

inline Vec4 Vec4::operator+(const Vec4& vector3) const
{
	return Vec4{ _m[0] + vector3._m[0], _m[1] + vector3._m[1], _m[2] + vector3._m[2], _m[3] + vector3._m[3] };
}

inline Vec4 Vec4::operator-(const Vec4& vector3) const
{
	return Vec4{ _m[0] - vector3._m[0], _m[1] - vector3._m[1], _m[2] - vector3._m[2], _m[3] - vector3._m[3] };
}

inline Vec4& Vec4::operator*=(float value)
{
	_m[0] *= value; _m[1] *= value; _m[2] *= value; _m[3] *= value;
	return *this;
}

inline Vec4& Vec4::operator/=(float value)
{
	_m[0] /= value; _m[1] /= value; _m[2] /= value; _m[3] /= value;
	return *this;
}

inline Vec4& Vec4::operator+=(const Vec4& vector3)
{
	_m[0] += vector3._m[0]; _m[1] += vector3._m[1]; _m[2] += vector3._m[2]; _m[3] += vector3._m[3];
	return *this;
}

inline Vec4& Vec4::operator-=(const Vec4& vector3)
{
	_m[0] -= vector3._m[0]; _m[1] -= vector3._m[1]; _m[2] -= vector3._m[2]; _m[3] -= vector3._m[3];
	return *this;
}

inline Vec4 Vec4::operator*(const Mat44& mat) const
{
	return Vec4(*this) *= mat;
}

inline Vec4& Vec4::operator*=(const Mat44& mat)
{
	Vec4 column{};
	Vec4 self = *this;
	mat.ToColumn(column, 0); _m[0] = self.Dot(column);
	mat.ToColumn(column, 1); _m[1] = self.Dot(column);
	mat.ToColumn(column, 2); _m[2] = self.Dot(column);
	mat.ToColumn(column, 3); _m[3] = self.Dot(column);
	return *this;
}

inline float Vec4::Dot(const Vec4& vector4) const
{
	return (_m[0] * vector4._m[0]) + (_m[1] * vector4._m[1]) + (_m[2] * vector4._m[2]) + (_m[3] * vector4._m[3]);
}

inline void Vec4::ToVec3(Vec3& output) const
{
	output._m[0] = _m[0];
	output._m[1] = _m[1];
	output._m[2] = _m[2];
}

inline void Mat44::ToColumn(Vec4& output, int column) const
{
	output._m[0] = _m[0][column];
	output._m[1] = _m[1][column];
	output._m[2] = _m[2][column];
	output._m[3] = _m[3][column];
}

inline void Mat44::ToRow(Vec4& output, int row) const
{
	output._m[0] = _m[row][0];
	output._m[1] = _m[row][1];
	output._m[2] = _m[row][2];
	output._m[3] = _m[row][3];
}

inline Mat44& Mat44::operator*=(const Mat44& mat)
{
	*this = operator*(mat);
	return *this;
}

inline Mat44 Mat44::operator*(const Mat44& mat) const
{
	Mat44 output;
	Vec4 row, column;
	//TODO: direct multiplication to reduce CPU overhead
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			ToRow(row, r); mat.ToColumn(column, c);
			output._m[r][c] = row.Dot(column);
		}
	}
	return output;
}

inline Mat44 Mat44::ProjectionLH(float width, float height, float fov, float zNear, float zFar)
{
	float aspectRatio = height / width;
	float invHalfTanTheta = 1.0f / tanf(fov * 3.14159f * 0.5f / 180.0f);
	float zScale = zFar / (zFar - zNear);

	Mat44 projectionMat{};
	projectionMat._m[0][0] = aspectRatio * invHalfTanTheta;		//projectionMat._m[0][1] = 0;				projectionMat._m[0][2] = 0;					projectionMat._m[0][3] = 0;
	/*	projectionMat._m[1][0] = 0;		*/							projectionMat._m[1][1] = invHalfTanTheta;	//projectionMat._m[1][2] = 0;				projectionMat._m[1][3] = 0;
	/*	projectionMat._m[2][0] = 0;									projectionMat._m[2][1] = 0;		*/			projectionMat._m[2][2] = zScale;			projectionMat._m[2][3] = 1;
	/*  projectionMat._m[3][0] = 0;									projectionMat._m[3][1] = 0;*/				projectionMat._m[3][2] = -zScale * zNear;	//projectionMat._m[3][3] = 0;
	return projectionMat;
}

inline Mat44 Mat44::ViewLH(const Vec3& eye, const Vec3& target, const Vec3& up)
{
	auto forward = (target - eye); forward.Normalize();
	auto right = up.Cross(forward); right.Normalize();
	auto upNew = forward.Cross(right); upNew.Normalize();

	Mat44 viewMat{};
	viewMat._m[0][0] = right.x;				viewMat._m[0][1] = upNew.x;			viewMat._m[0][2] = forward.x;//			viewMat._m[0][3] = 0;
	viewMat._m[1][0] = right.y;				viewMat._m[1][1] = upNew.y;			viewMat._m[1][2] = forward.y;//			viewMat._m[1][3] = 0;
	viewMat._m[2][0] = right.z;				viewMat._m[2][1] = upNew.z;			viewMat._m[2][2] = forward.z;//			viewMat._m[2][3] = 0;
	viewMat._m[3][0] = -right.Dot(eye);		viewMat._m[3][1] = -upNew.Dot(eye);	viewMat._m[3][2] = -forward.Dot(eye);	viewMat._m[3][3] = 1;
	return viewMat;
}

inline Mat44 Mat44::RotationX(float value)
{
	Mat44 rotMat{};
	rotMat._m[0][0] = 1;				/*rotMat._m[0][1] = 0;			rotMat._m[0][2] = 0;			rotMat._m[0][3] = 0;*/
	/*	rotMat._m[1][0] = 0;*/				rotMat._m[1][1] = cosf(value);	rotMat._m[1][2] = sinf(value);	//rotMat._m[1][3] = 0;
	/*	rotMat._m[2][0] = 0;*/				rotMat._m[2][1] = -sinf(value);	rotMat._m[2][2] = cosf(value);	//rotMat._m[2][3] = 0;
	/*	rotMat._m[3][0] = 0;				rotMat._m[3][1] = 0;			rotMat._m[3][2] = 0;*/			rotMat._m[3][3] = 1;
	return rotMat;
}

inline Mat44 Mat44::RotationY(float value)
{
	Mat44 rotMat{};
	rotMat._m[0][0] = cosf(value);/*	rotMat._m[0][1] = 0;*/			rotMat._m[0][2] = -sinf(value);	//rotMat._m[0][3] = 0;
	/*	rotMat._m[1][0] = 0;*/				rotMat._m[1][1] = 1;			//rotMat._m[1][2] = 0;			rotMat._m[1][3] = 0;
	rotMat._m[2][0] = sinf(value);		/*rotMat._m[2][1] = 0;*/		rotMat._m[2][2] = cosf(value);	//rotMat._m[2][3] = 0;
	/*	rotMat._m[3][0] = 0;*/				/*rotMat._m[3][1] = 0;			rotMat._m[3][2] = 0;*/			rotMat._m[3][3] = 1;
	return rotMat;
}

inline Mat44 Mat44::RotationZ(float value)
{
	Mat44 rotMat{};
	rotMat._m[0][0] = cosf(value);		rotMat._m[0][1] = sinf(value);	/*rotMat._m[0][2] = 0;	rotMat._m[0][3] = 0;*/
	rotMat._m[1][0] = -sin(value);		rotMat._m[1][1] = cos(value);	/*rotMat._m[1][2] = 0;	rotMat._m[1][3] = 0;*/
	/*	rotMat._m[2][0] = 0;				rotMat._m[2][1] = 0;*/			rotMat._m[2][2] = 1;	//rotMat._m[2][3] = 0;
	/*	rotMat._m[3][0] = 0;				rotMat._m[3][1] = 0;			rotMat._m[3][2] = 0;*/	rotMat._m[3][3] = 1;
	return rotMat;
}



