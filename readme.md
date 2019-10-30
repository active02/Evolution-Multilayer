Multilayer 반사 Spectrum의 계산은 Transfer Matrix Method 라는 Analytic 한 방법으로 가능하지만,
원하는 Spectrum을 갖는 Multilayer 구조를 찾는 일은 그리 간단한 일이 아닙니다.

본 파이썬 드는 이러한 Multilayer 구조의 최적화를 유전알고리즘을 이용하여 구현하기 위한 목적으로
2019년 2월에 작성하였던 프로그램입니다.

내용은 다음과 같습니다.
1. Multilayer의 각 층의 n, k, t 값을 주어진 범위에서 무작위로 생성하고, 생성된 구조변수에 해당하는 스펙트럼을 tmm 라이브러리를 이용하여 계산합니다.
2. 1의 과정으로 생성된 여러개의 데이터 중, 목표로하는 스펙트럼과 거리가 먼 것들을 제거한 후, 나머지 데이터들에 대해 crossover 및 mutation을 통해 차세대 구조변수들을 생성합니다. 여기에는 deap 라이브러리를 사용하였습니다.
3. 1-2 과정을 충분히 반복한 후 이를 시각화 합니다. 결과는 아래 그림처럼 나옵니다.






![Figure_1](https://github.com/active02/Evolution-Multilayer/blob/master/result.png)




그림설명

( 세대에 따른 구조 변화 : n 값)  ( 세대에 따른 구조 변화 : k 값)

( 최종 구조 : n 값)             ( 최종 구조 : k 값) ( x 축 : 깊이 (nm) )

(결과 스펙트럼)                 (세대에 따른 스펙트럼)

(세대에 따른 Target 스펙트럼에 대한 MAE 변화)
