# 蒙哥马利约分算法(Montgomery reduction algorithm, REDC)

蒙哥马利算法是进行快速连续模乘 `a * b % n` 的算法. 在这里, 我将详细解释算法的工作原理, 给出数学证明, 并提供示例代码作为演示. 蒙哥马利约分算法被广泛应用在密码学算法中, 例如 RSA 和 BN128.

## 求模运算

模运算即求余数, 例如 `10 % 3 = 1`. 模运算存在一些等式:

```text
(a + b) % n = (a % n + b % n) % n
(a - b) % n = (a % n - b % n) % n
(a * b) % n = (a % n * b % n) % n
```

## 逆元

若 `(a * b) % n = 1`, 则称 b 为 a 的模 n 逆元(inverse).

## 蒙哥马利约分

设 `N > 1`, 选一个 `R > N` 且 `gcd(R, N) = 1`. 则必然存在 `R⁻¹`, `N⁻¹` 满足:

```text
R * R⁻¹ % N = 1
R * R⁻¹ - N * N⁻¹ = 1
-N * N⁻¹ % R = 1
(R - N) * N⁻¹ % R = 1
0 < R⁻¹ < N
0 < N⁻¹ < R
```

定义蒙哥马利约分(Montgomery Reduction): `REDC(T) = TR⁻¹ % N`, 其中 `0 <= T < NR`, 则其通用计算过程伪代码如下:

```text
m = (T * N⁻¹) % R
t = (T + m * N) / R
if (N <= t):
    t -= N
```

注意到如果 R 是 2 的幂, 则除法和求模都可以使用移位完成.

## 蒙哥马利剩余

定义蒙哥马利剩余(Montgomery Residue) `X' = X * R % N`

```text
X' = X * R % N
   = X * R * R * R⁻¹ % N
   = REDC(X * (R * R % N))
```

```text
c  = a * b % N

a' = a * R % N
b' = b * R % N
c' = (a' * b') * R⁻¹ % N
   = REDC(a' * b')
c  = c' * R⁻¹ % n
   = REDC(c')
   = REDC(REDC(a' * b'))
```

其中 `R * R % N` 为常数, 在程序代码中可以进行提前计算.

## 代码实现

```py
import math

N = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47
R = 2 ** 256
assert R > N
assert math.gcd(N, R) == 1

R_INVERSE = 0x2e67157159e5c639cf63e9cfb74492d9eb2022850278edf8ed84884a014afa37
N_INVERSE = 0xf57a22b791888c6bd8afcbd01833da809ede7d651eca6ac987d20782e4866389
assert R * R_INVERSE % N == 1
assert R * R_INVERSE - N * N_INVERSE == 1
assert N_INVERSE * -N % R == 1
assert (R - N) * N_INVERSE % R == 1
assert 0 < R_INVERSE and R_INVERSE < N
assert 0 < N_INVERSE and N_INVERSE < R

R_POW2 = R * R % N


def redc(T):
    assert 0 <= T < R * N
    m = (T * N_INVERSE) % R
    t = (T + m * N) // R
    if t >= N:
        return t - N
    return t


def conv(x):
    return redc(x * R_POW2)


a = 0x1c658e925dbddaf46b81a8d835df5359f708114df717931be998b96a7fa69a18
b = 0x2f682d1f7dda8678b0d017978b3067b74807a5d49d2a41739659c6600a8bf018
assert (a * b) % N == redc(redc(conv(a) * conv(b)))


def fuzz():
    import random
    for _ in range(100000):
        a = random.randint(0, (2**256) - 1)
        b = random.randint(0, (2**256) - 1)
        assert (a * b) % N == redc(redc(conv(a) * conv(b)))
```
