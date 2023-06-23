from svg_turtle import SvgTurtle
from turtle import *


def koch(t,a,order):
    
    if order > 0:
        
        for alpha in [90,-90,-90,0,90,90,-90,0]:
            koch(t,a/4,order-1)
            t.left(alpha)
    else:
        t.forward(a)


def draw(t):

    size = 200
    t.penup()

    t.right(135)
    t.forward(size*1.4*0.5)
    t.right(-135)
    
    t.pendown()

    t.left(90)
    
    for i in range(4):
        koch(t,size,3)
        t.right(90)

def write_file(draw_func, filename, width, height):
    t = SvgTurtle(width, height)
    draw_func(t)
    t.save_as(filename)


def main():
    write_file(draw, 'example.svg', 800, 800)
    print('Done.')


if __name__ == '__main__':

    screen = Screen()
    screen.setup(500, 500)

    main()
    done()
