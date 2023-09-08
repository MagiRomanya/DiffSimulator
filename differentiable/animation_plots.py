#!/usr/bin/env python3

import matplotlib.pyplot as plt


class AnimatedPlot:
    """
    A straight forward way to animate plots using matplotlib.

    To use this class just start appending data to the plot with
    the get_data method.
    """

    def __init__(self, limdata=100):
        self.x_data = []
        self.y_data = []
        self.iteration = 0
        # limdata is an integer which tells the maximum number of stored points
        self.limdata = limdata
        self.fig, self.ax = plt.subplots()
        self.ax.set_title("Gradient decent convergence.")
        self.ax.set_xlabel("Iterations")
        self.ax.set_ylabel("Loss value")
        self.ax.grid()

    def get_data(self, value: float):
        """Append a new y value to the plot."""
        self.y_data.append(value)
        self.x_data.append(self.iteration)
        if len(self.y_data) > self.limdata:
            self.y_data = self.y_data[-self.limdata:]
            self.x_data = self.x_data[-self.limdata:]
        if (self.iteration == 0):
            self.__init_plot(value)

        self.iteration += 1

        self.line.set_data(self.x_data, self.y_data)
        self.ax.autoscale_view(True, True)
        self.ax.relim()
        self.fig.canvas.draw()
        plt.pause(0.00000000000001)

    def __init_plot(self, value):
        self.line, = self.ax.plot(value, 0, "x--")
