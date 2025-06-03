# -----------------------------------------------------------------------------
# 版权所有 (c) 2025 保留所有权利。
#
# 本文件隶属于 Poromechanics Solver 项目，主要开发者为：
#   - 马鹏飞：mapengfei@mail.nwpu.edu.cn
#   - 王璇：wangxuan2022@mail.nwpu.edu.cn
#
# 本软件仅供内部使用和学术研究之用。未经明确许可，严禁重新分发、修改或用于商业用途。
# 详细授权条款请参阅：https://www.pengfeima.cn/license-strict/
# -----------------------------------------------------------------------------

import numpy as np

# 输入参数
dt = 1e-3
R_AV = 0.15
R_MV = 0.375
R_sys = 27.0
C_AOR = 0.06
C_LA = 0.015

# TODO: 输入左心室的流量
P_AOR = np.loadtxt('p_ao.txt')
p_la = np.loadtxt('p_la.txt')
p_lv = np.loadtxt('p_lv.txt')
time_points = np.loadtxt('time.txt')
V_LV = np.loadtxt('v_lv.txt')
v_la = np.loadtxt('v_la.txt')

Q_LV = 0.0
Q_MV_n = 10.0
Q_sys_n = 20.0
Q_AV_n = 30.0
P_LA_n = 10.0
P_AOR_n = 100.0
P_LV_n = 100.0







def advance(Q_LV, Q_MV_n, Q_AV_n, Q_sys_n, P_LA_n, P_LV_n, P_AOR_n, dt):
    matrix = np.array([
        [1,     -1,   0,             0,               0,  0],
        [1,      0,   0,             0,         -1/R_MV,  1/R_MV],
        [0,      1,   0,        1/R_AV,               0, -1/R_AV],
        [0,      0,   1,      -1/R_sys,         1/R_sys,  0],
        [1,      0,   0,      -1/R_sys, C_LA/dt+1/R_sys,  0],
        [0,     -1,   0,C_AOR/dt+R_sys,        -1/R_sys,  0],
    ])
    matrix_inverse = np.linalg.inv(matrix)
    print(np.linalg.cond(matrix))

    Q_MV = Q_MV_n
    Q_AV = Q_AV_n
    Q_sys = Q_sys_n
    P_LA = P_LA_n
    P_LV = P_LV_n
    P_AOR = P_AOR_n
    rhs = np.array([Q_LV, 0, 0, 0, C_LA/dt*P_LA_n, C_AOR/dt*P_AOR_n])
    
    result = np.dot(matrix_inverse, rhs)
    Q_MV, P_LV, Q_AV, P_LA, P_AOR, Q_sys = result.flatten()
    # print(f"Q_MV = {Q_MV}, P_LV = {P_LV}, Q_AV = {Q_AV}, P_LA = {P_LA}, P_AOR = {P_AOR}, Q_sys = {Q_sys}")
    # print(Q_sys*Q_sys + Q_MV*Q_MV + Q_AV*Q_AV + P_LA*P_LA + P_AOR*P_AOR + P_LV*P_LV)
    return Q_sys, Q_MV, Q_AV, P_LA, P_AOR, P_LV




# 初始化存储列表
time_points = []
P_LV_list = []
P_AOR_list = []
P_LA_list = []
Q_LV_list = []
Q_sys_list = []
Q_AOR_list = []

for i in range(10000):
    time = i * dt
    Q_LV = 1000*np.cos(time)  # 模拟一个周期内的流量变化
    # Q_LV = time  # 假设左心室流量为常数
    Q_sys, Q_MV, Q_AV, P_LA, P_AOR, P_LV = advance(Q_LV, Q_MV_n, Q_AV_n, Q_sys_n, P_LA_n, P_LV_n, P_AOR_n, dt)
    
    # 更新下一时刻的值
    Q_MV_n = Q_MV
    Q_AV_n = Q_AV
    Q_sys_n = Q_sys
    P_LA_n = P_LA
    P_LV_n = P_LV
    P_AOR_n = P_AOR
    
    # 存储数据
    time_points.append(time)
    P_LV_list.append(P_LV)
    P_AOR_list.append(P_AOR)
    P_LA_list.append(P_LA)
    Q_LV_list.append(Q_LV)
    Q_sys_list.append(Q_sys)
    Q_AOR_list.append(Q_AV)
    # print(f"Step {i+1}: Q_sys = {Q_sys}, Q_MV = {Q_MV}, Q_AV = {Q_AV}, P_LA = {P_LA}, P_AOR = {P_AOR}, P_LV = {P_LV}")


import matplotlib.pyplot as plt
# 绘制图形
plt.figure(figsize=(18, 12))

# 绘制左心室压力
plt.subplot(2, 3, 1)
plt.plot(time_points, P_LV_list, label='P_LV (Left Ventricular Pressure)')
plt.xlabel('Time (s)')
plt.ylabel('Pressure (mmHg)')
plt.title('Left Ventricular Pressure over Time')
plt.grid(True)
plt.legend()

# 绘制主动脉压力
plt.subplot(2, 3, 2)
plt.plot(time_points, P_AOR_list, 'r', label='P_AOR (Aortic Pressure)')
plt.xlabel('Time (s)')
plt.ylabel('Pressure (mmHg)')
plt.title('Aortic Pressure over Time')
plt.grid(True)
plt.legend()

# 绘制左心房压力
plt.subplot(2, 3, 3)
plt.plot(time_points, P_LA_list, 'g', label='P_LA (Left Atrial Pressure)')
plt.xlabel('Time (s)')
plt.ylabel('Pressure (mmHg)')
plt.title('Left Atrial Pressure over Time')
plt.grid(True)
plt.legend()

# 绘制左心室流量
plt.subplot(2, 3, 4)
plt.plot(time_points, Q_LV_list, 'm', label='Q_LV (Left Ventricular Flow)')
plt.xlabel('Time (s)')
plt.ylabel('Flow (L/s)')
plt.title('Left Ventricular Flow over Time')
plt.grid(True)
plt.legend()

# 绘制左心室流量
plt.subplot(2, 3, 5)
plt.plot(time_points, Q_sys_list, 'm', label='Q_sys (Systemic Flow)')
plt.xlabel('Time (s)')
plt.ylabel('Flow (L/s)')
plt.title('Systemic Flow over Time')
plt.grid(True)
plt.legend()

# 绘制左心室流量
plt.subplot(2, 3, 6)
plt.plot(time_points, Q_AOR_list, 'm', label='Q_AOR (Aortic Flow)')
plt.xlabel('Time (s)')
plt.ylabel('Flow (L/s)')
plt.title('Aortic Flow over Time')
plt.grid(True)
plt.legend()

plt.tight_layout()
# plt.show()
plt.savefig("LV_hyper_dynamic_results.png")