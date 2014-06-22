#include <SDL2/SDL.h>
#include <memory>
#include "Game.h"

using namespace std;

// create game obj
std::shared_ptr<Game> g_game(new Game);

int main(int argc, char* argv[])
{
    g_game->init("Chapter 1", 100, 100, 640, 480, false);

    while (g_game->running())
    {
        g_game->handleEvents();
//        g_game->update();
        g_game->render();
    }
    g_game->clean();

    return 0;
}
