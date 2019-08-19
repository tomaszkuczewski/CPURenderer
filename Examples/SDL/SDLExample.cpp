// CPURenderer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#ifdef _MSC_VER
#pragma comment(lib, "SDL2.lib")
#endif

#define SDL_MAIN_HANDLED
#include <SDL.h>
#include <iostream>

#include "CPURenderer.h"

Mesh gCube{};

void initCube(Mesh* cube);


int main()
{
	initCube(&gCube);

	auto renderer = CPURenderer{640, 480};
	auto projection = Mat44::ProjectionLH(renderer.GetWidthf(), renderer.GetHeightf(), 90.0f, 0.1f, 1000.0f); //Projection matrix
	auto view = Mat44::ViewLH({ 2, 0, -1 }, { 0, 0, 1 }, Vec3{ 0, 1, 0 });

	renderer.SetProjectionMatrix(&projection);
	renderer.SetViewMatrix(&view);

	//Set light the of objects
	renderer.SetLight(Vec3{ 2, 0, -1 });

	//Create SDL Variables
	SDL_Renderer* sdlRenderer;
	SDL_Window* sdlWindow;
	SDL_Event sdlEvent;

	//Set the DRAW PIXEL LAMBDA
	renderer.SetDrawPixel([&](uint32_t x, uint32_t y, uint32_t color)->void {
		SDL_SetRenderDrawColor(sdlRenderer, (color >> 24) & 0xFF , (color >> 16) & 0xFF, (color >> 8) & 0xFF, color & 0xFF);
		SDL_RenderDrawPoint(sdlRenderer, x, y);
	});

	auto yaw = 0.0f;
	Mat44 rotation{};
	//Set the UPDATE LAMBDA
	renderer.SetUpdate([&](float delta)->void {
		//Handle events
		if (SDL_PollEvent(&sdlEvent) && sdlEvent.type == SDL_QUIT)
		{
			renderer.Stop();
			return;
		}
		//Clear screen
		renderer.Clear();
		//Apply world matrix with rotation
		rotation = Mat44::RotationY(yaw += delta);
		renderer.SetWorldMatrix(&rotation);
		//Draw
		renderer.Draw(&gCube);
		//As the buffer is ready then we can show it
		SDL_RenderPresent(sdlRenderer);
	});

	//Set the DRAW PIXEL LAMBDA
	renderer.SetDrawPixel([&](uint32_t x, uint32_t y, uint32_t color)->void {
		SDL_SetRenderDrawColor(sdlRenderer, (color >> 24) & 0xFF, (color >> 16) & 0xFF, (color >> 8) & 0xFF, color & 0xFF);
		SDL_RenderDrawPoint(sdlRenderer, x, y);
	});

	//SET THE START LAMBDA
	renderer.Start([&](int width, int height)->bool {
		if (SDL_Init(SDL_INIT_VIDEO) < 0)
			return false;

		if (SDL_CreateWindowAndRenderer(width, height, 0, &sdlWindow, &sdlRenderer))
			return false;

		return true;
	});
}

void initCube(Mesh* cube)
{
	cube->GetVertices().insert(cube->GetVertices().begin(),
	{	
		-0.5f, -0.5f, 0.5f, 1.0f, 0.0f, 0.0f, -1.0f,
		-0.5f, -0.5f, -0.5f, 1.0f, 0.0f, 0.0f, -1.0f,
		0.5f, -0.5f, -0.5f, 1.0f, 0.0f, 0.0f, -1.0f, 
		0.5f, -0.5f, -0.5f, 1.0f, 0.0f, 0.0f, -1.0f,
		0.5f, -0.5f, 0.5f, 1.0f, 0.0f, 0.0f, -1.0f,
		-0.5f, -0.5f, 0.5f, 1.0f, 0.0f, 0.0f, -1.0f,
		-0.5f, 0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 1.0f,
		0.5f, 0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 1.0f,
		0.5f, 0.5f, -0.5f, 1.0f, 0.0f, 0.0f, 1.0f,
		0.5f, 0.5f, -0.5f, 1.0f, 0.0f, 0.0f, 1.0f, 
		-0.5f, 0.5f, -0.5f, 1.0f, 0.0f, 0.0f, 1.0f, 
		-0.5f, 0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 1.0f, 
		-0.5f, -0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 0.0f,
		0.5f, -0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 0.0f,
		0.5f, 0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 0.0f, 
		0.5f, 0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 0.0f, 
		-0.5f, 0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 0.0f, 
		-0.5f, -0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 0.0f, 

		0.5f, -0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 0.0f, 
		0.5f, -0.5f, -0.5f, 1.0f, 0.0f, 0.0f, 0.0f, 
		0.5f, 0.5f, -0.5f, 1.0f, 0.0f, 0.0f, 0.0f, 
		0.5f, 0.5f, -0.5f, 1.0f, 0.0f, 0.0f, 0.0f,
		0.5f, 0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 0.0f, 
		0.5f, -0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 0.0f, 
		0.5f, -0.5f, -0.5f, 1.0f, 0.0f, 0.0f, 0.0f, 
		-0.5f, -0.5f, -0.5f, 1.0f, 0.0f, 0.0f, 0.0f,
		-0.5f, 0.5f, -0.5f, 1.0f, 0.0f, 0.0f, 0.0f,
		-0.5f, 0.5f, -0.5f, 1.0f, 0.0f, 0.0f, 0.0f,
		0.5f, 0.5f, -0.5f, 1.0f, 0.0f, 0.0f, 0.0f,
		0.5f, -0.5f, -0.5f, 1.0f, 0.0f, 0.0f, 0.0f, 
		-0.5f, -0.5f, -0.5f, 1.0f, 0.0f, 0.0f, 0.0f,
		-0.5f, -0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 0.0f,
		-0.5f, 0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 0.0f,
		-0.5f, 0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 0.0f, 
		-0.5f, 0.5f, -0.5f, 1.0f, 0.0f, 0.0f, 0.0f,
		-0.5f, -0.5f, -0.5f, 1.0f, 0.0f, 0.0f, 0.0f
	});

	for (int i = 0; i < 36; i++)
		cube->GetIndices().push_back(i);
}
