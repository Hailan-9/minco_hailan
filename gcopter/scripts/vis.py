#!/usr/bin/env python3
# coding=utf-8
from cProfile import label
import rospy
import math
from matplotlib import pyplot as plt
from std_msgs.msg import Float64MultiArray
import yaml
import rospkg
class Vis:
    def __init__(self):

        self.bs_sub = rospy.Subscriber(
            "/visualizer/bs", Float64MultiArray, self.bs_callback)
        self.as_sub = rospy.Subscriber(
            "/visualizer/As", Float64MultiArray, self.as_callback)
        self.s_sub = rospy.Subscriber(
            "/visualizer/s", Float64MultiArray, self.s_callback)
        self.accx_sub = rospy.Subscriber(
            "/visualizer/accx", Float64MultiArray, self.accx_callback)
        self.accy_sub = rospy.Subscriber(
            "/visualizer/accy", Float64MultiArray, self.accy_callback)
        self.v_sub = rospy.Subscriber(
            "/visualizer/v", Float64MultiArray, self.v_callback)

        self.bs = []
        self.As = []
        self.s = []
        self.accx = []
        self.accy = []
        self.v = []

        self.rcv_bs = False
        self.rcv_as = False
        self.rcv_s = False
        self.rcv_accx = False
        self.rcv_accy = False
        self.rcv_v = False
        rospack = rospkg.RosPack()
        yamlFileDir = rospack.get_path('gcopter') + '/config/curve_gen.yaml'
        with open(yamlFileDir, 'r') as y:
            self.yaml = yaml.load(y)
            self.v_max = self.yaml['v_max']
            self.a_max = self.yaml['a_max']

        self.id = 1;

    def start(self):
        rospy.init_node("vis")
        rate = rospy.Rate(1)
        print("start")

        while not rospy.is_shutdown():

            if self.rcv_s and self.rcv_v and self.rcv_accx and self.accy:
                # figsize是画布尺寸，facecolor是背景颜色，num为图形编号或者名称
                plt.figure(num=self.id, figsize=(2, 3), facecolor='gray')
                self.id = self.id + 1
                l1 = plt.plot(self.s, self.v, color='blue',
                              linewidth=1.0, label='speed')
                l2 = plt.plot(self.s, self.accx, color='red',
                              linewidth=1.0, label='acc_x')
                l3 = plt.plot(self.s, self.accy, color='green',
                              linewidth=1.0, label='acc_y')
                l4 = plt.plot([self.s[0], self.s[-1]],
                              [self.v_max, self.v_max], color='lightsalmon', linewidth=1.0, linestyle='--', label='v_max')
                l5 = plt.plot([self.s[0], self.s[-1]],
                              [self.a_max, self.a_max], color='lightpink', linewidth=1.0, linestyle='--', label='a_max')
                l6 = plt.plot([self.s[0], self.s[-1]],
                              [-self.a_max, -self.a_max], color='lightpink', linewidth=1.0, linestyle='--', label='a_mim')
                plt.legend()
                plt.xlabel("s")
                plt.show()
                self.rcv_s = False
                self.rcv_accx = False
                self.rcv_accy = False
                self.bs = []
                self.s = []
                self.accx = []
                self.accy = []
                self.v = []

            self.rcv_bs = False
            self.rcv_as = False
            self.rcv_v = False

            rate.sleep()

    def bs_callback(self, msg):
        self.rcv_bs = True
        self.rcv_v = True
        n = msg.layout.dim[0].size
        for i in range(n):
            self.bs.append(msg.data[i])
            self.v.append(math.sqrt(msg.data[i]))

    def as_callback(self, msg):
        self.rcv_as = True
        n = msg.layout.dim[0].size
        for i in range(n):
            self.As .append(msg.data[i])

    def s_callback(self, msg):
        print("recive s")

        self.rcv_s = True
        n = msg.layout.dim[0].size
        for i in range(n):
            self.s.append(msg.data[i])

    def accx_callback(self, msg):
        self.rcv_accx = True
        n = msg.layout.dim[0].size
        for i in range(n):
            self.accx.append(msg.data[i])
        self.accx.append(0.0)

    def accy_callback(self, msg):
        self.rcv_accy = True
        n = msg.layout.dim[0].size
        for i in range(n):
            self.accy.append(msg.data[i])
        self.accy.append(0.0)

    def v_callback(self, msg):
        self.rcv_v = True
        n = msg.layout.dim[0].size
        for i in range(n):
            self.v.append(msg.data[i])


if __name__ == "__main__":
    vis = Vis()
    vis.start()
