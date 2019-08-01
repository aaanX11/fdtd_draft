import numpy as np
import matplotlib.pyplot as plt

tau = 1e-6
def f(t):
    return 300 + 300 * (1 - np.exp(-t / tau))
def f_pr(t):
    return 300 * np.exp(-t / tau) / tau
def f_prpr(t):
    return -300 * np.exp(-t / tau) / tau ** 2

def f1(t):
    return np.sin(np.pi*(t-1e-10)/1e-10)
def f1_pr(t):
    return 1e+10*np.pi*np.cos(np.pi*(t-1e-10)/1e-10)
def f1_prpr(t):
    return -1e+20*np.pi*np.pi*np.sin(np.pi*(t-1e-10)/1e-10)


def test1():
    tau = 1e-6
    def f(t):
        return 300 + 300*(1-np.exp(-t/tau))
    def f_pr(t):
        return 300*np.exp(-t/tau)/tau
    def f_prpr(t):
        return -300*np.exp(-t/tau)/tau**2

    t = 0

    T = [f(t)]

    s1 = f_pr(t)
    s_ = [s1]
    T_ = [s1]

    s2 = f_prpr(t)
    T__ = [s2]
    s__ = [s2]

    dt = 2e-14
    a_ = []
    for i in range(8):
        t += dt
        T.append(f(t))
        T_.append(f_pr(t))
        T__.append(f_prpr(t))

        a = (T[-1] - T[-2]) / dt ** 3 - s1 / dt ** 2 - 0.5 * s2 / dt
        a_.append(a)

        s1_ = 3*(T[-1] - T[-2]) / dt - 2*s1 - 0.5*s2*dt
        s2_ = 6*(T[-1] - T[-2]) / dt**2 - 6*s1/dt - 2*s2
        s1 = s1_
        s2 = s2_

        s_.append(s1)
        s__.append(s2)


    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    ax1.set_title('Real derivative')
    ax1.plot(T_)

    ax2.set_title('Num derivative')
    ax2.plot(s_)

    ax3.set_title('Real 2 derivative')
    ax3.plot(T__)
    print T__

    ax4.plot(s__)
    print s__

    print a_

    plt.show()
    raw_input()


def test2():

    t = 0
    dt = 2e-14

    T = [f1(-dt), f1(t)]

    s1 = f1_pr(t)
    s_ = [s1]
    T_ = [s1]

    s2 = f1_prpr(t)
    T__ = [s2]
    s__ = [s2]


    a_ = []
    for i in range(100):
        t += dt
        T.append(f1(t))
        T_.append(f1_pr(t))
        T__.append(f1_prpr(t))

        a = (T[-1] - T[-2]) / dt ** 3 - s1 / dt ** 2 - 0.5 * s2 / dt
        a_.append(a)

        s1 += s2*dt + 3*a*dt**2
        s2 += 6*a*dt

        if abs(s1) > 10*abs(T_[-1]) or abs(s2) > 10*abs(T__[-1]):
            print 'i = ', i
            break

        s_.append(s1)
        s__.append(s2)

    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    ax1.set_title('Real derivative')
    ax1.plot(T_)

    ax2.set_title('Num derivative')
    ax2.plot(s_)

    ax3.set_title('Real 2 derivative')
    ax3.plot(T__)
    print T__

    ax4.plot(s__)
    print s__

    print a_

    plt.show()
    raw_input()


def test3():

    t = 0
    dt = 2e-14

    T = [f1(-dt), f1(t)]

    s1 = f1_pr(t)
    s_ = [s1]
    T_ = [s1]

    s2 = f1_prpr(t)
    T__ = [s2]
    s__ = [s2]


    a_ = []
    for i in range(100):
        t += dt
        T.append(f1(t))
        T_.append(f1_pr(t))
        T__.append(f1_prpr(t))

        a = (T[-1] - T[-2]) / dt ** 3 - s1 / dt ** 2 - 0.5 * s2 / dt
        a_.append(a)

        #s1 += s2*dt + 3*a*dt**2
        #s2 += 6*a*dt

        s2 = (T[-1] -2*T[-2] + T[-3])/dt**2
        s1_ = 3*(T[-1] - T[-2]) / dt - 2*s1 - 0.5*s2*dt
        s1 = s1_

        if abs(s1) > 10*abs(T_[-1]) or abs(s2) > 10*abs(T__[-1]):
            print 'i = ', i
            break

        s_.append(s1)
        s__.append(s2)

    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    ax1.set_title('Real derivative')
    ax1.plot(T_)
    print 'real der'
    print T_

    ax2.set_title('Num derivative')
    ax2.plot(s_)
    print 'num der'
    print s_

    ax3.set_title('Real 2 derivative')
    ax3.plot(T__)
    print 'real 2der'
    print T__

    ax4.plot(s__)
    print 'num 2der'
    print s__

    plt.show()
    raw_input()


def test4():

    t = 0
    dt = 2e-14

    T = [f1(-dt), f1(t)]

    s1 = f1_pr(t)
    s_ = [s1]
    T_ = [s1]

    s2 = f1_prpr(t)
    T__ = [s2]
    s__ = []

    a_ = []
    for i in range(100):
        t += dt
        T.append(f1(t))
        T_.append(f1_pr(t))
        T__.append(f1_prpr(t))

        s2 = (T[-1] - 2*T[-2] + T[-3])/dt**2
        s1 = (T[-1] - T[-3]) / (2*dt)

        a = (T[-1] - T[-2]) / dt ** 3 - s1 / dt ** 2 - 0.5 * s2 / dt
        a_.append(a)

        if abs(s1) > 10*abs(T_[-1]) or abs(s2) > 10*abs(T__[-1]):
            print 'i = ', i
            break

        s_.append(s1)
        s__.append(s2)

    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    ax1.set_title('Real derivative')
    ax1.plot(T_)

    ax2.set_title('Num derivative')
    ax2.plot(s_)

    ax3.set_title('Real 2 derivative')
    ax3.plot(T__)
    print 'real der2'
    print T__

    ax4.set_title('Num 2 derivative')
    ax4.plot(s__)
    print 'num der2'
    print s__

    plt.show()
    raw_input()


def test5():
    t = 0
    dt = 2e-14

    T = [f1(-dt), f1(t)]

    s1 = f1_pr(t)
    s_ = [s1]
    T_ = [s1]

    s2 = f1_prpr(t)
    T__ = [s2]
    s__ = [s2]

    a_ = []
    for i in range(100):
        t += dt
        T.append(f1(t))
        T_.append(f1_pr(t))
        T__.append(f1_prpr(t))

        # s1 += s2*dt + 3*a*dt**2
        # s2 += 6*a*dt

        s2 = (T[-1] - 2 * T[-2] + T[-3]) / dt ** 2
        s1_ = 3 * (T[-1] - T[-2]) / dt - 2 * s1 - 0.5 * s2 * dt
        # s2_ = 6*(T[-1] - T[-2]) / dt**2 - 6*s1/dt - 2*s2
        s1 = s1_
        # s2 = s2_

        a = (T[-1] - T[-2]) / dt ** 3 - s1 / dt ** 2 - 0.5 * s2 / dt
        a_.append(a)

        if abs(s1) > 10 * abs(T_[-1]) or abs(s2) > 10 * abs(T__[-1]):
            print 'i = ', i
            break

        s_.append(s1)
        s__.append(s2)

    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    ax1.set_title('Real derivative')
    ax1.plot(T_)

    ax2.set_title('Num derivative')
    ax2.plot(s_)

    ax3.set_title('Real 2 derivative')
    ax3.plot(T__)
    print T__

    ax4.plot(s__)
    print s__

    print a_

    plt.show()
    raw_input()