import matplotlib.pyplot as plt

from turtle import *

screen = Screen()
screen.setup(500, 500)


size = 200

penup()

right(135)
forward(size*1.4*0.5)
right(-135)

pendown()

def koch(a, order):
    
    if order > 0:
        
        for t in [90,-90,-90,0,90,90,-90,0]:
            koch(a/4,order-1)
            left(t)
    else:
        forward(a)


left(90)
tracer(100)
hideturtle()

for i in range(4):
    koch(size,3)
    right(90)


update()

done()


