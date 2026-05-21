---
tags:
  - probabilidades
  - estadistica
  - transformaciones
  - casella-berger
fecha: 2026-05-13
---
[[00. Probabilidades I]]

> [!abstract] Material Complementario
> Desarrollo detallado de los ejemplos y teoremas del libro **Casella & Berger** sobre transformaciones bivariadas. Complementa los apuntes de la [[13. Transformaciones Bivariadas (Casella & Berger)|Semana 13]]. Utiliza conceptos de [[06. Esperanza, Varianza y Función Generadora de Momentos|función generadora de momentos]], [[10. Distribuciones Bivariadas - Dist. Conjunta, Marginales, Varianza y Correlación|distribuciones conjuntas y marginales]], [[11. Distribuciones Bivariadas - Dist. Condicionales e Independencia|independencia]] y [[09. Transformación de Variables Aleatorias Univariadas|transformaciones univariadas]]. Las distribuciones empleadas incluyen la [[08. Modelos de Probabilidad Continuos|Normal, Exponencial y Beta]] y la [[07. Modelos de Probabilidad Discretos|Poisson]].

---

**Ejemplo 4.2.13 (Fgm de una suma de variables normales)** A veces el Teorema 4.2.12 puede usarse para derivar fácilmente la distribución de $Z$ a partir del conocimiento de las distribuciones de $X$ e $Y$. Por ejemplo, sea $X \sim n(\mu, \sigma^2)$ e $Y \sim n(\gamma, \tau^2)$ variables aleatorias normales independientes. Del Ejercicio 2.33, las fgm de $X$ e $Y$ son

$$M_X(t) = \exp(\mu t + \sigma^2t^2/2) \quad \text{y} \quad M_Y(t) = \exp(\gamma t + \tau^2t^2/2).$$

Así, del Teorema 4.2.12, la fgm de $Z = X + Y$ es

$$M_Z(t) = M_X(t)M_Y(t) = \exp \left( (\mu + \gamma)t + (\sigma^2 + \tau^2)t^2/2 \right).$$

Esta es la fgm de una variable aleatoria normal con media $\mu + \gamma$ y varianza $\sigma^2 + \tau^2$. Este resultado es lo suficientemente importante como para enunciarse como teorema. $\square$

> [!theorem] Teorema 4.2.14 _Sea $X \sim n(\mu, \sigma^2)$ e $Y \sim n(\gamma, \tau^2)$ variables aleatorias normales independientes. Entonces la variable aleatoria $Z = X + Y$ tiene distribución $n(\mu + \gamma, \sigma^2 + \tau^2)$._

Si $f(x, y)$ es la fdp conjunta del vector aleatorio continuo $(X, Y)$, la condición (4.2.1) puede no cumplirse en un conjunto $A$ de valores $(x, y)$ para el cual $\int_A \int dx , dy = 0$. En tal caso, $X$ e $Y$ se siguen denominando variables aleatorias independientes. Esto refleja el hecho de que dos fdp que difieren únicamente en un conjunto como $A$ definen la misma distribución de probabilidad para $(X, Y)$. Para ver esto, supóngase que $f(x, y)$ y $f^_(x, y)$ son dos fdp que coinciden en todos lados excepto en un conjunto $A$ para el cual $\int_A \int dx , dy = 0$. Sea $(X, Y)$ con fdp $f(x, y)$, sea $(X^_, Y^_)$ con fdp $f^_(x, y)$, y sea $B$ cualquier subconjunto de $\mathbb{R}^2$. Entonces

$$ \begin{aligned} P((X, Y) \in B) &= \int \int_B f(x, y) , dx , dy \ &= \int \int_{B \cap A^c} f(x, y) , dx , dy \ &= \int \int_{B \cap A^c} f^_(x, y) , dx , dy \ &= \int \int_B f^_(x, y) , dx , dy = P((X^_, Y^_) \in B). \end{aligned} $$

Así, $(X,Y)$ y $(X^_,Y^_)$ tienen la misma distribución de probabilidad. Por ejemplo, $f(x,y) = e^{-x-y}$, $x > 0$ e $y > 0$, es una fdp para dos variables aleatorias exponenciales independientes y satisface (4.2.1). Pero $f^_(x,y)$, que es igual a $f(x,y)$ excepto que $f^_(x,y) = 0$ si $x = y$, también es la fdp de dos variables aleatorias exponenciales independientes, aunque (4.2.1) no se cumpla en el conjunto $A = {(x,x) : x > 0}$.

---

## 4.3 Transformaciones Bivariadas

En la Sección 2.1 se discutieron métodos para encontrar la distribución de una función de una variable aleatoria. En esta sección extendemos esas ideas al caso de vectores aleatorios bivariados.

Sea $(X,Y)$ un vector aleatorio bivariado con distribución de probabilidad conocida. Consideremos ahora un nuevo vector aleatorio bivariado $(U,V)$ definido por $U = g_1(X,Y)$ y $V = g_2(X,Y)$, donde $g_1(x,y)$ y $g_2(x,y)$ son funciones especificadas. Si $B$ es cualquier subconjunto de $\mathbb{R}^2$, entonces $(U,V) \in B$ si y solo si $(X,Y) \in A$, donde $A = {(x,y) : (g_1(x,y), g_2(x,y)) \in B}$. Así, $P((U,V) \in B) = P((X,Y) \in A)$, y la distribución de probabilidad de $(U,V)$ queda completamente determinada por la distribución de probabilidad de $(X,Y)$.

Si $(X,Y)$ es un vector aleatorio bivariado discreto, entonces existe solo un conjunto numerable de valores para los cuales la fmp conjunta de $(X,Y)$ es positiva. Llamemos a este conjunto $\mathcal{A}$. Definamos el conjunto $B = {(u,v) : u = g_1(x,y) \text{ y } v = g_2(x,y) \text{ para algún } (x,y) \in \mathcal{A}}$. Entonces $B$ es el conjunto numerable de valores posibles para el vector aleatorio discreto $(U,V)$. Y si, para cualquier $(u,v) \in B$, se define $A_{uv}$ como ${(x,y) \in \mathcal{A} : g_1(x,y) = u \text{ y } g_2(x,y) = v}$, entonces la fmp conjunta de $(U,V)$, $f_{U,V}(u,v)$, puede calcularse a partir de la fmp conjunta de $(X,Y)$ mediante

$$(4.3.1) \quad f_{U,V}(u,v) = P(U = u, V = v) = P((X,Y) \in A_{uv}) = \sum_{(x,y) \in A_{uv}} f_{X,Y}(x,y).$$

**Ejemplo 4.3.1 (Distribución de la suma de variables Poisson)** Sean $X$ e $Y$ variables aleatorias de Poisson independientes con parámetros $\theta$ y $\lambda$, respectivamente. Así, la fmp conjunta de $(X,Y)$ es

$$f_{X,Y}(x,y) = \frac{\theta^x e^{-\theta}}{x!} \frac{\lambda^y e^{-\lambda}}{y!}, \quad x = 0, 1, 2, \dots, , y = 0, 1, 2, \dots$$

El conjunto $\mathcal{A}$ es ${(x, y) : x = 0, 1, 2, \dots \text{ y } y = 0, 1, 2, \dots}$. Definamos ahora $U = X + Y$ y $V = Y$, es decir, $g_1(x,y) = x + y$ y $g_2(x,y) = y$. Describamos el conjunto $\mathcal{B}$, el conjunto de posibles valores $(u,v)$. Los valores posibles para $v$ son los enteros no negativos, ya que $v = y$. Para un valor fijo de $V = v$, se tiene que $u = x + y = x + v$ debe ser un entero mayor o igual a $v$, dado que $x$ es un entero no negativo. El conjunto de todos los posibles valores $(u,v)$ es entonces $\mathcal{B} = {(u, v) : v = 0, 1, 2, \dots \text{ y } u = v, v + 1, v + 2, \dots}$. Para cualquier $(u, v) \in \mathcal{B}$, el único valor $(x,y)$ que satisface $x + y = u$ e $y = v$ es $x = u - v$ e $y = v$. Por lo tanto, en este ejemplo $A_{uv}$ consiste siempre en el único punto $(u - v, v)$. De (4.3.1) obtenemos entonces la fmp conjunta de $(U,V)$ como

$$f_{U,V}(u,v) = f_{X,Y}(u-v, v) = \frac{\theta^{u-v}e^{-\theta}}{(u-v)!} \frac{\lambda^v e^{-\lambda}}{v!}, \quad \begin{array}{l} v = 0, 1, 2, \dots, \ u = v, v+1, v+2, \dots \end{array}$$

En este ejemplo es interesante calcular la fmp marginal de $U$. Para cualquier entero no negativo fijo $u$, se tiene que $f_{U,V}(u,v) > 0$ únicamente para $v = 0, 1, \dots, u$. Esto nos da el conjunto de valores de $v$ sobre los cuales sumar para obtener la fmp marginal de $U$:

$$f_U(u) = \sum_{v=0}^u \frac{\theta^{u-v}e^{-\theta}}{(u-v)!} \frac{\lambda^v e^{-\lambda}}{v!} = e^{-(\theta+\lambda)} \sum_{v=0}^u \frac{\theta^{u-v} \lambda^v}{(u-v)! v!}, \quad u = 0, 1, 2, \dots$$

Esto puede simplificarse notando que, si multiplicamos y dividimos cada término por $u!$, podemos aplicar el Teorema del Binomio para obtener

$$f_U(u) = \frac{e^{-(\theta+\lambda)}}{u!} \sum_{v=0}^u \binom{u}{v} \lambda^v \theta^{u-v} = \frac{e^{-(\theta+\lambda)}}{u!} (\theta+\lambda)^u, \quad u = 0, 1, 2, \dots$$

Esta es la fmp de una variable aleatoria de Poisson con parámetro $\theta + \lambda$. Este resultado es lo suficientemente significativo como para enunciarse como teorema.

> [!theorem] Teorema 4.3.2 _Si $X \sim \text{Poisson}(\theta)$ e $Y \sim \text{Poisson}(\lambda)$ son independientes, entonces $X + Y \sim \text{Poisson}(\theta + \lambda)$._

Si $(X, Y)$ es un vector aleatorio continuo con fdp conjunta $f_{X,Y}(x, y)$, entonces la fdp conjunta de $(U, V)$ puede expresarse en términos de $f_{X,Y}(x, y)$ de manera análoga a (2.1.8). Como antes, $A = {(x, y) : f_{X,Y}(x, y) > 0}$ y $B = {(u, v) : u = g_1(x, y) \text{ y } v = g_2(x, y) \text{ para algún } (x, y) \in A}$. La fdp conjunta $f_{U,V}(u, v)$ será positiva en el conjunto $B$. Para la versión más simple de este resultado, suponemos que la transformación $u = g_1(x, y)$ y $v = g_2(x, y)$ define una transformación biyectiva de $A$ sobre $B$. La transformación es sobreyectiva por la definición de $B$. Estamos suponiendo que para cada $(u, v) \in B$ existe un único $(x, y) \in A$ tal que $(u, v) = (g_1(x, y), g_2(x, y))$. Para tal transformación biyectiva, podemos resolver las ecuaciones $u = g_1(x, y)$ y $v = g_2(x, y)$ para expresar $x$ e $y$ en términos de $u$ y $v$. Denotaremos esta transformación inversa por $x = h_1(u, v)$ e $y = h_2(u, v)$. El papel que desempeña la derivada en el caso univariado es ahora desempeñado por una cantidad denominada el _Jacobiano de la transformación_. Esta función de $(u, v)$, denotada por $J$, es el determinante de una _matriz_ de derivadas parciales, definida por

$$ J = \begin{vmatrix} \dfrac{\partial x}{\partial u} & \dfrac{\partial x}{\partial v} \[8pt] \dfrac{\partial y}{\partial u} & \dfrac{\partial y}{\partial v} \end{vmatrix} = \frac{\partial x}{\partial u} \frac{\partial y}{\partial v} - \frac{\partial y}{\partial u} \frac{\partial x}{\partial v}, $$

donde

$$\frac{\partial x}{\partial u} = \frac{\partial h_1(u, v)}{\partial u}, \quad \frac{\partial x}{\partial v} = \frac{\partial h_1(u, v)}{\partial v}, \quad \frac{\partial y}{\partial u} = \frac{\partial h_2(u, v)}{\partial u}, \quad \text{y} \quad \frac{\partial y}{\partial v} = \frac{\partial h_2(u, v)}{\partial v}.$$

Suponemos que $J$ no es idénticamente $0$ en $B$. Entonces la fdp conjunta de $(U, V)$ es $0$ fuera del conjunto $B$, y sobre el conjunto $B$ está dada por

$$(4.3.2) \qquad\qquad\qquad\qquad f_{U,V}(u, v) = f_{X,Y}(h_1(u, v), h_2(u, v))|J|,$$

donde $|J|$ es el valor absoluto de $J$. Al aplicar (4.3.2), a veces resulta igual de difícil determinar el conjunto $B$ y verificar que la transformación es biyectiva que sustituir directamente en la fórmula (4.3.2). Obsérvense estos aspectos en los siguientes ejemplos.

**Ejemplo 4.3.3 (Distribución del producto de variables beta)** Sea $X \sim \text{beta}(\alpha, \beta)$ e $Y \sim \text{beta}(\alpha + \beta, \gamma)$ variables aleatorias independientes. La fdp conjunta de $(X, Y)$ es

$$f_{X,Y}(x, y) = \frac{\Gamma(\alpha + \beta)}{\Gamma(\alpha)\Gamma(\beta)} x^{\alpha - 1}(1 - x)^{\beta - 1} \frac{\Gamma(\alpha + \beta + \gamma)}{\Gamma(\alpha + \beta)\Gamma(\gamma)} y^{\alpha + \beta - 1}(1 - y)^{\gamma - 1},$$

$$0 < x < 1, \quad 0 < y < 1.$$

Consideremos la transformación $U = XY$ y $V = X$. El conjunto de posibles valores de $V$ es $0 < v < 1$, ya que $V = X$. Para un valor fijo de $V = v$, $U$ debe estar entre $0$ y $v$, puesto que $X = V = v$ e $Y$ está entre $0$ y $1$. Así, esta transformación mapea el conjunto $A$ sobre el conjunto $B = {(u, v) : 0 < u < v < 1}$. Para cualquier $(u, v) \in B$, las ecuaciones $u = xy$ y $v = x$ pueden resolverse de forma única para obtener $x = h_1(u, v) = v$ e $y = h_2(u, v) = u/v$. Nótese que, considerada como transformación definida sobre todo $\mathbb{R}^2$, esta transformación no es inyectiva: cualquier punto $(0, y)$ se mapea al punto $(0, 0)$. Pero como función definida únicamente sobre $A$, sí es una transformación biyectiva sobre $B$. El Jacobiano está dado por

$$ J = \begin{vmatrix} \dfrac{\partial x}{\partial u} & \dfrac{\partial x}{\partial v} \[8pt] \dfrac{\partial y}{\partial u} & \dfrac{\partial y}{\partial v} \end{vmatrix} = \begin{vmatrix} 0 & 1 \[6pt] \dfrac{1}{v} & -\dfrac{u}{v^2} \end{vmatrix} = -\frac{1}{v}. $$

Así, de (4.3.2) obtenemos la fdp conjunta como

$$(4.3.3) \quad f_{U,V}(u, v) = \frac{\Gamma(\alpha + \beta + \gamma)}{\Gamma(\alpha)\Gamma(\beta)\Gamma(\gamma)} v^{\alpha - 1}(1 - v)^{\beta - 1} \left(\frac{u}{v}\right)^{\alpha + \beta - 1} \left(1 - \frac{u}{v}\right)^{\gamma - 1} \frac{1}{v},$$

$$\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad 0 < u < v < 1.$$

La distribución marginal de $V = X$ es, naturalmente, una distribución $\text{beta}(\alpha, \beta)$. Pero la distribución de $U$ también es una distribución beta:

$$f_U(u) = \int_u^1 f_{U,V}(u, v) , dv$$

$$= \frac{\Gamma(\alpha + \beta + \gamma)}{\Gamma(\alpha)\Gamma(\beta)\Gamma(\gamma)} u^{\alpha - 1} \int_u^1 \left( \frac{u}{v} - u \right)^{\beta - 1} \left( 1 - \frac{u}{v} \right)^{\gamma - 1} \left( \frac{u}{v^2} \right) , dv.$$

Se usó la expresión (4.3.3) pero algunos términos han sido reordenados. Ahora realizamos el cambio de variable univariado $y = (u/v - u)/(1 - u)$, de modo que $dy = -u/[v^2(1 - u)]dv$, para obtener

$$ \begin{aligned} f_U(u) &= \frac{\Gamma(\alpha + \beta + \gamma)}{\Gamma(\alpha)\Gamma(\beta)\Gamma(\gamma)} u^{\alpha - 1} (1 - u)^{\beta + \gamma - 1} \int_0^1 y^{\beta - 1} (1 - y)^{\gamma - 1} , dy \ &= \frac{\Gamma(\alpha + \beta + \gamma)}{\Gamma(\alpha)\Gamma(\beta)\Gamma(\gamma)} u^{\alpha - 1} (1 - u)^{\beta + \gamma - 1} \frac{\Gamma(\beta)\Gamma(\gamma)}{\Gamma(\beta + \gamma)} \ &= \frac{\Gamma(\alpha + \beta + \gamma)}{\Gamma(\alpha)\Gamma(\beta + \gamma)} u^{\alpha - 1} (1 - u)^{\beta + \gamma - 1}, \quad 0 < u < 1. \end{aligned} $$

Para obtener la segunda identidad reconocemos el integrando como el núcleo de una fdp beta y usamos (3.3.17). Así, la distribución marginal de $U$ es $\text{beta}(\alpha, \beta + \gamma)$. $\square$

**Ejemplo 4.3.4 (Suma y diferencia de variables normales)** Sean $X$ e $Y$ variables aleatorias normales estándar independientes. Consideremos la transformación $U = X + Y$ y $V = X - Y$. En la notación usada anteriormente, $U = g_1(X, Y)$ donde $g_1(x, y) = x + y$ y $V = g_2(X, Y)$ donde $g_2(x, y) = x - y$. La fdp conjunta de $X$ e $Y$ es $f_{X,Y}(x, y) = (2\pi)^{-1} \exp(-x^2/2) \exp(-y^2/2)$, con $-\infty < x < \infty$ y $-\infty < y < \infty$. Así el conjunto $\mathcal{A} = \mathbb{R}^2$. Para determinar el conjunto $\mathcal{B}$ sobre el cual $f_{U,V}(u,v)$ es positiva, debemos determinar todos los valores que

$$(4.3.4) \qquad\qquad\qquad\qquad u = x + y \quad \text{y} \quad v = x - y$$

toman cuando $(x, y)$ recorre el conjunto $\mathcal{A} = \mathbb{R}^2$. Pero podemos fijar $u$ y $v$ como cualquier número y resolver de forma única las ecuaciones (4.3.4) para $x$ e $y$, obteniendo

$$(4.3.5) \qquad\qquad x = h_1(u,v) = \frac{u + v}{2} \quad \text{y} \quad y = h_2(u,v) = \frac{u - v}{2}.$$

Esto demuestra dos cosas. Para cualquier $(u, v) \in \mathbb{R}^2$ existe un $(x, y) \in \mathcal{A}$ (definido por (4.3.5)) tal que $u = x + y$ y $v = x - y$, por lo que $\mathcal{B} = \mathbb{R}^2$. Además, dado que la solución (4.3.5) es única, esto también muestra que la transformación considerada es biyectiva. De (4.3.5) las derivadas parciales de $x$ e $y$ son fáciles de calcular:

$$ J = \begin{vmatrix} \dfrac{\partial x}{\partial u} & \dfrac{\partial x}{\partial v} \[8pt] \dfrac{\partial y}{\partial u} & \dfrac{\partial y}{\partial v} \end{vmatrix} = \begin{vmatrix} \dfrac{1}{2} & \dfrac{1}{2} \[6pt] \dfrac{1}{2} & -\dfrac{1}{2} \end{vmatrix} = -\frac{1}{2}. $$

Sustituyendo las expresiones (4.3.5) para $x$ e $y$ en $f_{X,Y}(x, y)$ y usando $|J| = \frac{1}{2}$, obtenemos la fdp conjunta de $(U, V)$ a partir de (4.3.2) como

$$f_{U,V}(u, v) = f_{X,Y}(h_1(u, v), h_2(u, v))|J| = \frac{1}{2\pi} e^{-((u+v)/2)^2/2} e^{-((u-v)/2)^2/2} \frac{1}{2}$$

para $-\infty < u < \infty$ y $-\infty < v < \infty$. Al desarrollar los cuadrados en las exponenciales, los términos con $uv$ se cancelan, y tras simplificar y reordenar se obtiene

$$f_{U,V}(u, v) = \left( \frac{1}{\sqrt{2\pi}\sqrt{2}} e^{-u^2/4} \right) \left( \frac{1}{\sqrt{2\pi}\sqrt{2}} e^{-v^2/4} \right).$$

La fdp conjunta se factorizó en una función de $u$ y una función de $v$. Por el Lema 4.2.7, $U$ y $V$ son independientes. Por el Teorema 4.2.14, la distribución marginal de $U = X+Y$ es $n(0, 2)$. Análogamente, el Teorema 4.2.12 puede usarse para mostrar que la distribución marginal de $V$ también es $n(0, 2)$. Este hecho importante —que sumas y diferencias de variables normales independientes son también variables normales independientes— es válido independientemente de las medias de $X$ e $Y$, siempre que $\text{Var}(X) = \text{Var}(Y)$. Este resultado se deja como Ejercicio 4.27. Los Teoremas 4.2.12 y 4.2.14 nos dan las distribuciones marginales de $U$ y $V$, pero el análisis más detallado aquí realizado es necesario para determinar que $U$ y $V$ son independientes. $\square$

En el Ejemplo 4.3.4 encontramos que $U$ y $V$ son variables aleatorias independientes. Existe una situación mucho más simple, pero muy importante, en la que nuevas variables $U$ y $V$ definidas en términos de las variables originales $X$ e $Y$ resultan ser independientes. El Teorema 4.3.5 describe esta situación.