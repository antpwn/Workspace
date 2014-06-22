#include <SDL2/SDL.h>

SDL_Window* g_pWindow = 0;
SDL_Renderer* g_pRenderer = 0;

int main(int argc, char* args[])
{
    // initialize SDL and check
    if (SDL_Init(SDL_INIT_EVERYTHING) >= 0)
    {
        // succeded. lets create our window
        g_pWindow = SDL_CreateWindow("Chapter 1: Setting up SDL",
                    SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                    640, 480,
                    SDL_WINDOW_SHOWN);

        // check for window creation
        if (g_pWindow != 0)
        {
            g_pRenderer = SDL_CreateRenderer(g_pWindow, -1, 0);
        }
    }
    else
    {
        return 1;
    }

    //start drawing
    // set to black, function accepts R,G,B,Alpha
    SDL_SetRenderDrawColor(g_pRenderer, 255, 0, 0, 255);

    // clear window to black
    SDL_RenderClear(g_pRenderer);

    // show window
    SDL_RenderPresent(g_pRenderer);

    // delay before exit
    SDL_Delay(5000);

    // clean up SDL
    SDL_Quit();

    // exit
    return 0;
}
