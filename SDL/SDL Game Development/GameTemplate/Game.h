#ifndef __Game__
#define __Game__

class Game
{
public:
    Game() {}
    ~Game() {}

    // set running var to true
    bool init(const char* title,
              int xpos, int ypos,
              int width, int height,
              bool fullscreen);

    void render();
    void update();
    void handleEvents();
    void clean();

    // function to access the private running var
    bool running() { return m_bRunning; }

private:
    bool m_bRunning;

    SDL_Window* m_pWindow;
    SDL_Renderer* m_pRenderer;
};

#endif // __Game__
