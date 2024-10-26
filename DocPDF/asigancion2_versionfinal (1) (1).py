import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy import constants

# Definimos la ecuación diferencial en coordenadas polares
def V(r, lamda):
    return (k * np.exp(-lamda * r)) / r

# Definir la función para r_min usando `r_min_alternativo`
def r_min_alternativo(r, b, e, lamda):
    return (r**2) - (b**2) - ((k * r * np.exp(-lamda * r)) / e)

# Definir el integrando de la expresión para theta_max
def integrando(r, b, e, lamda):
    # Potencial de Yukawa
    potential = V(r, lamda)
    expr = 1 - (potential / e) - (b**2 / r**2)
    if expr <= 0:
        return np.inf
    return b / (r**2 * np.sqrt(expr))

# Definir la función de ángulo de dispersión
def Angulo_dispersion(theta_max):
    return np.pi - 2 * theta_max

def seccion_eficaz(Parametros_b, Chis, db_dchis):
   exp = (Parametros_b/np.sin(Chis))*np.abs(db_dchis)
   return exp

def Valor_critico(L, m, k):
   exp = (-1)*(1.1905*L**2)/(m*k)
   return exp


# Definimos las condiciones iniciales 
q1 = constants.e * 79  # Masa que esta en el foco (Nucleo de oro)
q2 = constants.e * 2   # Particula Alpha
m = 6.64424e-27       # Masa de la particula alpha
k = (q1*q2)/(4*np.pi*constants.epsilon_0)  # Constante del potencial
print(k)
b = 0.1  # Parámetro de impacto   # A tener en cuenta 1
lamda = 1 # Parametro de apantallamiento de Yukawa
L = 1 #Momento Angular
r0 = 10 # Posición radial inicial
dr0 = -5 # Velocidad radial inicial
theta0 = 3*np.pi/2  # Ángulo inicial

def V(r, lamda):
    return (k * np.exp((-1)*r*lamda)) / r  #Potencial Repulsivo 
    
# Definimos el potencial efectivo
def V_ef(r, lamda):
    return (V(r, lamda)) + ((L**2)/(2*m*(r**2))) 

#Grafico del potencial ####################################################################3

x_vals = np.linspace(1e-40, 1, 100000)
f_vals = V_ef(x_vals, lamda)

plt.ylim([max(f_vals.min(), 1e-5), f_vals.max()]) 
plt.plot(x_vals, f_vals, label=r'$Potencial$')
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Potencial de Efectivo de Yukawa')
plt.yscale('log')
plt.legend()
plt.show()


Punto_critico = Valor_critico(L, 1, 1)
print("Punto critico", Punto_critico)
#Grafica de los potenciales para distintos valores de lambda, tomando el vlaor critico
m = 1
k = 1
#Ya se tiene definida una malla en lineas anterior, se reutiliza
x_vals = np.linspace(0.01, 10, 1000)
Potencial_evaluado_1 = V_ef(x_vals, Punto_critico+(Punto_critico/100))
Potencial_evaluado_2 = V_ef(x_vals, Punto_critico)
Potencial_evaluado_3 = V_ef(x_vals, -Punto_critico+(Punto_critico/100))

plt.ylim([0, 10])
plt.xlim([0, 5])
# Primer potencial
plt.plot(x_vals, Potencial_evaluado_1, label=f'Potencial 1 lambda igual a {Punto_critico+(Punto_critico/100)}', color='blue')
# Segundo potencial con color diferente
plt.plot(x_vals, Potencial_evaluado_2, label=f'Potencial 2 lambda igual a {Punto_critico}', color='green')
# Tercer potencial con otro color y estilo de línea
plt.plot(x_vals, Potencial_evaluado_3, label=f'Potencial 3 lambda igual a {-Punto_critico+(Punto_critico/100)}', color='red')

# Añadir etiquetas, leyenda y detalles de la gráfica
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Potenciales Efectivos de Yukawa')
plt.legend()

# Mostrar la gráfica
plt.show()

######################################################################

#Grafica energias vs angulo de dispersion

#Para un parametro alfa dado, varias la energia


#######################################################################
print("1. Ejemplo")
print("2. Simulación de ángulos de dispersión y determinación del ángulo solido")
print("3. Simulación para n particulas y gráfica de Energía vs Ángulo de dispersión para distintos λ")
print("4. Simulación de ángulos de dispersión y determinación del ángulo solido")
print("5. Simulación y gráfica de Sección Eficaz vs Ángulo de Dispersión")
menu = int(input("Ingrese 1, 2, 3, 4 o 5: "))

if menu == 1:
 q1 = constants.e * 79  # Masa que esta en el foco (Nucleo de oro)
 q2 = constants.e * 2   # Particula Alpha
 m = 6.64424e-27       # Masa de la particula alpha
 k = (q1*q2)/(4*np.pi*constants.epsilon_0)  # Constante del potencial
 print(k)
 b = 0.0001  # Parámetro de impacto   # A tener en cuenta 1
 lamda = 1 # Parametro de apantallamiento de Yukawa
 L = 1 #Momento Angular
 r0 = 10 # Posición radial inicial
 dr0 = -5 # Velocidad radial inicial
 theta0 = 3*np.pi/2  # Ángulo inicial
 # Tiempo de simulación

 #Al final, se termina utilizando la forma "alternativa", dado a problemas con el metodo numerico que se utiliza al calcular la energia para cada r con valores de las constantes muy pequeñas
 e = (0.5)*m*dr0**2
 x_initial_guess = [0, 10] #Este intervalo debe ser modificado respecto a la observación de la gráfica anterior
 raices_1 = [fsolve(r_min_alternativo, x_initial_guess, args=(b, e, lamda))[0], fsolve(r_min_alternativo, x_initial_guess, args=(b, e, lamda))[1]]

 for i in range(0,2):
    if raices_1[i] > 0:
       distancia_min = raices_1[i]
 print("La distancia mínima alternativa es igual a : ", distancia_min)

# Realizamos la integración numérica desde r_min hasta "infinito", osea, 1000
 theta_max, error = quad(integrando, distancia_min, 1000, args=(b, e, lamda)) #Se utiliza la función 'quad' de integracion por cuadratura de gauss para integrar la función "integrand" que se creo en el apartado anterior, es decir, hallamos tetha maximo 

# Mostramos el resultado
 print(f"El valor de theta_max es: {theta_max:.4f}")
 print(f"Error estimado de la integración: {error:.4e}")




#Calculo del angulo de dispersión Chi
 Angulo_de_dispersion = Angulo_dispersion(theta_max)

#Lo anterior solo fue un ejemplo..

 print("EL angulo de dispersión es: ", Angulo_de_dispersion)


if menu == 2:
 print("La idea principal de esta 'simulación', es calcular el ángulo de dispersión Chi para un conjunto de particulas que tengan un parametro de impacto comprendido dentro de un intervalo definido, para posteriormente /n Encontral el ángulo solido al que fueron dispersadas y encontrar la sección eficaz. ")
 #Parametros Iniciales
 numero_de_particulas = 100
 lista_distancias_minimas = np.zeros(numero_de_particulas)
 lista_angulos_maximos = np.zeros(numero_de_particulas)
 lista_chi_angulos = np.zeros(numero_de_particulas) #Angulos de dispersion
 k = 1 # Valor arbitrario
 m = 1 #masa
 Parametros_impacto = np.random.uniform(0.1, 1, numero_de_particulas) #generamos un valor b "aleatorio" para cada particula
 lamda = 1 # Parametro de apantallamiento de Yukawa
 L = 1 #Momento Angular

 r0 = 10 # Posición radial inicial
 dr0 = -5 # Velocidad radial inicial
 theta0 = 3*np.pi/2  # Ángulo inicial
 Area = 5
 x_initial_guess = [0, 10] #Intervalo en donde se buscan las raices
 #Se calcula la energia del sistema 
 e = (0.5)*m*dr0**2
 print("La energía total del sistema es aproximadamente :", e)
 
 for i in range(numero_de_particulas):
    b = Parametros_impacto[i]
    raices = [fsolve(r_min_alternativo, x_initial_guess, args=(b, e, lamda))[0], fsolve(r_min_alternativo, x_initial_guess, args=(b, e, lamda))[1]]
    for x in range(0,2):
      if raices[x] > 0:
       lista_distancias_minimas[i] = raices[x]
    lista_angulos_maximos[i], error = quad(integrando, lista_distancias_minimas[i], 1000, args=(b, e, lamda))
    lista_chi_angulos[i] = np.pi -  2*lista_angulos_maximos[i]

 lista_ordenada_de_angulos_chi = np.sort(lista_chi_angulos)

 print("todo bien", lista_chi_angulos)

 print(lista_ordenada_de_angulos_chi[0], lista_ordenada_de_angulos_chi[-1]) #Muestra el angulo menor y mayor

 #Calculo del angulo solido

 angulo_solido = 2*np.pi*((-np.cos(lista_ordenada_de_angulos_chi[-1])) - (-np.cos(lista_ordenada_de_angulos_chi[0])))

if menu == 3:
 numero_de_particulas = 100
 #Constantes
 q1 = constants.e * 79  # Masa que esta en el foco (Nucleo de oro)
 q2 = constants.e * 2   # Particula Alpha
 m = 6.64424e-27       # Masa de la particula alpha
 k = (q1*q2)/(4*np.pi*constants.epsilon_0)  # Constante del potencial
 print(k)
 k= 1
 m=1
 b = 0.01  # Parámetro de impacto   # A tener en cuenta 1
 lamda = 1 # Parametro de apantallamiento de Yukawa
 valores_lambda = [0.5, 1, 5, 10]  # Ajustar según lo que quieras
 colores = ['blue', 'green', 'red', 'purple']  # Colores para cada lambda
 L = 1 #Momento Angular
 r0 = 10 # Posición radial inicial
 dr0 = np.linspace(1, 1000, numero_de_particulas)# Velocidad radial inicial
 theta0 = 3*np.pi/2  # Ángulo inicial
 x_initial_guess = [0, 10] #Intervalo en donde se buscan las raices
 marcadores = ['o', 's', '^', 'D']  # Diferentes marcadores para cada lambda

 plt.figure(figsize=(12, 8))
 for j, lamda in enumerate(valores_lambda):
    lista_distancias_minimas = np.zeros(numero_de_particulas)
    lista_angulos_maximos = np.zeros(numero_de_particulas)
    lista_chi_angulos = np.zeros(numero_de_particulas)  # Angulos de dispersion
    lista_energias = []

    for i in range(numero_de_particulas):
        e = (0.5) * m * (dr0[i]) ** 2  # Energía cinética
        lista_energias.append(e)
        
        # Encontrar las raíces de la ecuación para obtener r_min
        raices = fsolve(r_min_alternativo, x_initial_guess, args=(b, e, lamda))
        raices_validas = [r for r in raices if r > 0]  # Solo usar raíces positivas

        if len(raices_validas) > 0:
            lista_distancias_minimas[i] = raices_validas[0]  # Asignar la raíz válida
        
        # Calcular ángulo máximo mediante integración
        lista_angulos_maximos[i], error = quad(integrando, lista_distancias_minimas[i], 1000, args=(b, e, lamda))
        lista_chi_angulos[i] = np.pi - 2 * lista_angulos_maximos[i]  # Ángulo de dispersión

    # Graficar para este valor de lambda con líneas y marcadores diferentes
    plt.plot(lista_chi_angulos, lista_energias, color=colores[j], marker=marcadores[j], label=f'λ = {lamda}', linestyle='-', markersize=5)

# Configurar los ejes y la gráfica
 plt.xlabel('Ángulo de dispersión (radianes)', fontsize=14)
 plt.ylabel('Energía (Joules)', fontsize=14)
 plt.title('Energía vs Ángulo de dispersión para distintos λ', fontsize=16)
 plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# Configuración de la escala
 plt.yscale('log')  # Escala logarítmica para energía

# Añadir divisiones a los ejes
 plt.xticks(fontsize=12)
 plt.yticks(fontsize=12)
 plt.minorticks_on()  # Activar las divisiones menores
 plt.grid(True, which='minor', linestyle=':', linewidth=0.5)

# Añadir leyenda
 plt.legend(fontsize=12)

# Mostrar la gráfica
 plt.show()

##################################################################
if menu == 4:
 #Grafica de Ángulo de Dispersión vs Parámetro de Impacto para Diferentes λ 
 numero_de_particulas = 100
 m = 1 # Masa de la partícula alpha
 k = 1 # Constante del potencial de Yukawa
 L = 1 # Momento angular
 e = 1 #1 julio, 1 MeV

# Parámetros de impacto y de apantallamiento
 b_values = np.linspace(0.01, 5, numero_de_particulas)
 lambda_values = [0, 1, 2, 3, 5]
 

 plt.figure(figsize=(10, 6))

 for lamda in lambda_values:
    angulos_dispersion = []
    for b in b_values:
        e = 1

        # Calcular r_min usando fsolve con r_min_alternativo
        rmin_solution = fsolve(r_min_alternativo, 1, args=(b, e, lamda))
        rmin = rmin_solution[0] if rmin_solution[0] > 0 else np.nan
        
        # Si rmin es válido, calcular el ángulo de dispersión
        if not np.isnan(rmin):
            try:
                # Calcular theta_max
                theta_max, _ = quad(integrando, rmin, np.inf, args=(b, e, lamda))
                chi = Angulo_dispersion(theta_max)
                angulos_dispersion.append(chi)
            except Exception as error:
                print(f"Error en la integración para b={b}, λ={lamda}: {error}")
                angulos_dispersion.append(np.nan)
        else:
            angulos_dispersion.append(np.nan)

    # Graficar ángulo de dispersión vs parámetro de impacto para cada lambda
    plt.plot(b_values, angulos_dispersion, label=f'λ = {lamda}')

# Configuración del gráfico
 plt.xlabel('Parámetro de Impacto (b)', fontsize=14)
 plt.ylabel('Ángulo de Dispersión (rad)', fontsize=14)
 plt.title('Ángulo de Dispersión vs Parámetro de Impacto para Diferentes λ', fontsize=16)
 plt.grid(True, which='both', linestyle='--', linewidth=0.5)
 plt.legend()
 plt.show()

 
# Definir constantes y parámetros iniciales
numero_de_particulas = 100  # Número de partículas
q1 = constants.e * 79       # Carga del núcleo (Núcleo de oro)
q2 = constants.e * 2        # Carga de la partícula (Partícula Alpha)
m = 6.64424e-27             # Masa de la partícula alpha
k = (q1 * q2) / (4 * np.pi * constants.epsilon_0)  # Constante del potencial de Yukawa
b = np.linspace(0.01, 1, numero_de_particulas)  # Parámetro de impacto
lamda = 1                # Parámetro de apantallamiento de Yukawa
L = 1                    # Momento angular
dr0 = np.linspace(1, 1000, numero_de_particulas)  # Velocidad radial inicial

# Listas para almacenar los resultados
lista_angulos_dispersion = np.zeros(numero_de_particulas)
lista_secciones_eficaces = np.zeros(numero_de_particulas)


if menu == 5:
 numero_de_particulas = 100  # Número de partículas
 q1 = constants.e * 79       # Carga del núcleo (Núcleo de oro)
 q2 = constants.e * 2        # Carga de la partícula (Partícula Alpha)
 m = 6.64424e-27             # Masa de la partícula alpha
 k = (q1 * q2) / (4 * np.pi * constants.epsilon_0)  # Constante del potencial de Yukawa
 b_values = np.linspace(0.01, 1, numero_de_particulas)  # Parámetro de impacto
 lamda = 1                   # Parámetro de apantallamiento de Yukawa
 L = 1                       # Momento angular
 dr0 = np.linspace(1, 1000, numero_de_particulas)  # Velocidad radial inicial


# Listas para almacenar los resultados
 lista_angulos_dispersion = np.zeros(numero_de_particulas)
 lista_secciones_eficaces = np.zeros(numero_de_particulas)

# Cálculo del ángulo de dispersión y sección eficaz para cada partícula
 for i in range(numero_de_particulas):
    e = 0.5 * m * dr0[i]**2  # Energía inicial

    # Resolver r_min usando fsolve con r_min_alternativo
    rmin_solution = fsolve(r_min_alternativo, 1, args=(b_values[i], e, lamda))
    rmin = rmin_solution[0] if rmin_solution[0] > 0 else np.nan

    # Si rmin es válido, calcular el ángulo de dispersión
    if not np.isnan(rmin):
        theta_max, _ = quad(integrando, rmin, np.inf, args=(b_values[i], e, lamda))
        chi = Angulo_dispersion(theta_max)
        lista_angulos_dispersion[i] = chi
    else:
        lista_angulos_dispersion[i] = np.nan

# Calcular la derivada db/dchi
 db_dchis = np.gradient(b_values, lista_angulos_dispersion)

# Cálculo de la sección eficaz
 lista_secciones_eficaces = seccion_eficaz(b_values, lista_angulos_dispersion, db_dchis)

# Graficar la sección eficaz vs ángulo de dispersión
 plt.figure(figsize=(10, 6))
 plt.plot(lista_angulos_dispersion, lista_secciones_eficaces, 'o-', color='blue', alpha=0.7, label=r'$\sigma(\chi)$')

# Configurar los ejes y el gráfico
 plt.xlabel('Ángulo de Dispersión (rad)', fontsize=14)
 plt.ylabel('Sección Eficaz', fontsize=14)
 plt.title('Sección Eficaz vs Ángulo de Dispersión', fontsize=16)
 plt.yscale('log')  # Escala logarítmica en el eje y si es necesario
 plt.grid(True, which='both', linestyle='--', linewidth=0.5)
 plt.legend()
 plt.show()