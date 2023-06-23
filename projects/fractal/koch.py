import matplotlib.pyplot as plt

from turtle import *

screen = Screen()
screen.setup(500, 500)


size = 200

penup()

right(135)
forward(size*0.5)
right(-135)

pendown()

def koch(a, order):
    
    if order > 0:
        
        for t in [60, -120, 60, 0]:
            koch(a/3,order-1)
            left(t)
    else:
        forward(a)

left(60)

tracer(100)
hideturtle()

for i in range(3):
    koch(size,4)
    right(120)

update()
done()
