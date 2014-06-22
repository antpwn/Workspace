#include <SDL2/SDL.h>

// SDLs
SDL_Window* g_pWindow = 0;
SDL_Renderer* g_pRenderer = 0;

// loop bool
bool g_bRunning = false;

bool init(const char* title,
          int xpos, int ypos,
          int height, int width,
          int flags)
{
    // initialize SDL and check
    if (SDL_Init(SDL_INIT_EVERYTHING) >= 0)
    {
        // succeded. lets create our window
        g_pWindow = SDL_CreateWindow(title,
                    xpos, ypos,
                    height, width,
                    flags);

        // check for window creation
        if (g_pWindow != 0)
        {
            g_pRenderer = SDL_CreateRenderer(g_pWindow, -1, 0);
        }
    }
    else
    {
        return false;
    }
    return true;
}

void render()
{
    // set to black, function accepts R,G,B,Alpha
    SDL_SetRenderDrawColor(g_pRenderer, 255, 0, 0, 255);

    // clear window to black
    SDL_RenderClear(g_pRenderer);

    // show window
    SDL_RenderPresent(g_pRenderer);
}

int main(int argc, char* args[])
{
    if (init("Chapter 1: Setting up SDL",
             SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
             640, 480,
             SDL_WINDOW_SHOWN))
    {
        g_bRunning = true;
    }
    else
    {
        return 1;
    }

    //start drawing
    while (g_bRunning)
    {
        render();
    }

    // clean up SDL
    SDL_Quit();

    // exit
    return 0;
}


