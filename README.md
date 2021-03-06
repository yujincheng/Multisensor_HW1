# Multi-sensor Homework 1

本readme文档即本次作业的实验报告，要求实现教材第11章中2个算例中的一个，要求分别采用Dolev算法、快速收敛算法、传感器融合算法及混合算法等4种算法。

源码在sensor.m文件中，报告参考本文件，或者 ./doc/README.pdf

## 拜占庭将军问题的基本定理

### 离散的拜占庭将军定理
有m个将军，他们从指挥官那里得到命令会转发给别人进行相互验证。但是，他们之间存在叛徒，并且可以相互通信，叛徒的指令是不可预测的。那么如果想让所有忠诚将军都按照指挥官的计划行动，则，最多存在k个叛徒，m>3k

### 连续的拜占庭将军问题

有m个传感器，他们测量同一个标量，并且把自己的测量值转发给其他人。其中有些传感器是错误的，他们转发的值不可预期。那么，同样的，如果想共同获得一个比较精确的估计，则，最多只能k个传感器坏了，m>3k。

### 任务

本次作业需要验证的算法就是，在m=5，k=1的时候，用几个算法来计算4个正常传感器应该输出的值。

## Dovlev 算法

Dolev 算法的思路实在是非常简单，总而言之概括成为一句话，去掉最高分，去掉最低分。如果坏的传感器刚好被去掉了，那么结果是可信的。如果坏的传感器没被去掉，那么说明坏的传感器的值刚好还落到了正常传感器测量范围内。

dolve算法最简单，不存在迭代。需要排序和计算平均值。

在每个机器上，只需要统计最大值和最小值。同时计算平均。所以单个机器上的计算复杂度为。

计算量瓶颈在于计算t个最小值，本算法复杂度为 $min(O(nt),O(nlog(n))$

## 快速收敛算法

快速收敛法让人难以理解。

听名字是一个迭代算法，但实际上看算法描述就是把超出一定范围的数用范围内一个数取代。

## 传感器融合算法

传感器融合算法相比Dolev更近一步了。因为数据都是传感器多次测量得到的。单个传感器测量本身就会存在很大的不确定度。所以会统计数据范围内有多少传感器能够覆盖。如果一个范围内，有至少 m-k 个传感器能够测量，则说明该范围是合理的。

传感器融合算法难点在于统计传感器区间个数困难，并且统计区间个数成为算法的瓶颈。每个机器需要对 2N个数（包括N个机器的最大最小值）进行排序。

计算量瓶颈在这里与这个排序，算法复杂度为 $O(2nlog(2n)) = O(nlog(n))$

## 混合算法

相比于传感器融合算法，找可信传感器区间的思路。混合算法的传感器更进一步。想通过统计的区间计算出来一个确切值来。思想就是各个可信区间的值进行加权。

我实现的方案很简单，用区间的平均数来代表这个区间，用可信区间内的传感器测量数目来作为该区间的权重。利用该代表元和权重进行算术加权求和。

混合算法仅需要在传感器融合算法计算可信区间的平均数，增加的计算量可以忽略。

所以算法复杂度仍为 $O(nlog(n))$。
