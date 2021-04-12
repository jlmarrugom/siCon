import streamlit as st
from scipy.integrate import odeint
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#st.set_page_config(layout="wide")
st.set_page_config(page_title='SiConApp',page_icon=':100:')
st.title('Simulador de Pandemias :earth_americas:')
@st.cache
def deriv(y, t, N, beta, alfaa, alfag, alfah, alfai, deltaa, deltai, gammaa, gammai, gammah, gammag, sigmah, sigmag, pe, pa, pi, ph, pg, w, mu, v, lamda, q):
    #Ecuaciones del cambio para solucionar
    S, E, A, I, H, G, R, V, FV, D = y
    dSdt = (-beta(t) * S * (pe*E+pa*A+pi*I+ph*H+pg*G) / N) -(v(t)*w*S) +lamda -(mu*S) #v(t)*S vacunados en tiempo t, solo vacunación sobre susceptibles es efectiva
    dEdt = (beta(t) * (S+FV) * (pe*E+pa*A+pi*I+ph*H+pg*G) / N) - ((deltaa+deltai+mu) * E)
    dAdt = deltaa*E - (gammaa+alfaa+mu)*A
    dIdt = deltai * E - (sigmah +sigmag+gammai+alfai+mu)*I
    dHdt = sigmah*I - (gammah+alfah+mu)*H#*mu # el mu hace que crezca mucho
    dGdt = sigmag*I - (gammag+alfag+mu)*G
    dRdt = gammaa*A + gammai * I + gammah*H + gammag*G - mu*R
    dVdt = v(t)*q*w*S- mu*V #las vacunaciones tienen una componente temporal
    dFVdt= v(t)*(1-q)*w*S - (beta(t)/N)*FV*(pe*E+pa*A+pi*I+ph*H+pg*G) - mu*FV #comp temporal
    N = S+E+A+I+H+G+R+V+FV
    dDdt = mu * N
    return dSdt, dEdt, dAdt, dIdt, dHdt, dGdt, dRdt, dVdt, dFVdt, dDdt

#La sidebar NO se actualiza cada vez que recarga el archivo!!!
N = st.sidebar.number_input('Población Total',min_value=1_000,max_value=100_000_000,value=1_000_000,step=1_000)

R_0 = st.sidebar.number_input('Tasa de Reproducción Base (R0)',min_value=0.0,max_value=10.0,value=2.5,step=0.5)

tap = st.sidebar.slider('Porcentaje Uso del tapabocas',0,100,0,1)
tap = tap*0.35/100

hig = st.sidebar.slider('Porcentaje de Higiene (lavado de manos)',0,100,0,1)
hig = hig*0.25/100

lamda = st.sidebar.number_input('Tasa de Natalidad',value=0.01)
alfag = st.sidebar.number_input('Tasa de Mortalidad inducida por enfermedad para pacientes de UCI',value=0.01)
sigmag = st.sidebar.number_input('Tasa de Ingreso a UCI',value=0.01)

#beta = (1-tap)*(1-hig)*R_0 * gamma   # R_0 = beta / gamma, so beta = R_0 * gamma


if st.sidebar.checkbox('¿Hay Encerramiento?'):

    L = st.sidebar.number_input('Comienzo del Encerramiento',min_value=1,value=200,step=1)
    dur = st.sidebar.slider('Duración',min_value=1,max_value=60,value=20,step=1)
else: 
    L = 1000
    dur=0
#Para el encerramiento:
def R__0(t):
    return R_0 if ((t < L)|(t>(L+dur))) else 0.9
def beta(t):
    return (1-tap)*(1-hig)*R__0(t)

if st.sidebar.checkbox('¿Hay vacunación?'):
    vac = st.sidebar.number_input('Tasa de Vacunación',value=5*10e-3)
    dia = st.sidebar.number_input('Día de Vacunación',1,None,10,1)
    q = st.sidebar.number_input('Eficacia',value=1)
    w = st.sidebar.number_input('Adhesión a la campaña',value=1)
else:
    vac=0
    dia=0
    q=0
    w=0
def v(t):
    return vac if (t<dia+1+vac/100)&(t>=dia) else 0.0   #Tasa de vacunación

t_max = st.sidebar.number_input('Tiempo Máximo de Simulación',min_value=100,value=400,step=10)

# st.write('Parametros: ','\u03C3:',str(sigma),'\u03B3:',str(gamma),'\u03B2:',str(beta(0)),'\u03BE:',str(xi),'\u03BC',str(mu),'\u03BD',str(v(0)))
S0, E0, A0, I0, H0, G0, R0, V0, FV0, D0 = N-1 , 1, 0, 0, 0, 0, 0, 0, 0, 0  # initial conditions: one exposed
#S, E, A, I, H, G, R, V, FV, D
#A los susceptibles les quitamos los vacunados, y a los recuperados se los sumamos (En las ecuaciones)
t = np.linspace(0, t_max, t_max)
#Es necesario que sean t_max pasos para que el df sea acorde

y0 = S0, E0, A0, I0, H0, G0, R0, V0, FV0, D0 # Initial conditions vector
alfaa, alfah, alfai = 0, 10e-4, 5*10e-4 #tasa de mortalidad
deltaa, deltai = 0.874, 1/2 #Tasa de evolución asin, infec
gammaa, gammai, gammah, gammag = 1/14,1/3,1/10,0.05 #Tasa de Recuperación Asintomaticos, infectados, Hospital,UCI
sigmah = 0.12 #Tasa de Hospitalización
pe, pa, pi, ph, pg = 0.4, 1/3, 1.0, 0.01, 0.01 #Infecciosidad de personas Expuestas
mu = 3.91*10e-5 #Tasa de Mortalidad natural

# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, alfaa, alfag, alfah, alfai, deltaa, deltai, gammaa, gammai, gammah, gammag, sigmah, sigmag, pe, pa, pi, ph, pg, w, mu, v, lamda, q))
S, E, A, I, H, G, R, V, FV, D = ret.T

population = pd.DataFrame(columns=['Susceptibles','Expuestos','Infectados','Recuperados','Fallecidos','Total'])
population['Susceptibles'] = S
population['Expuestos'] = E
population['Asintomáticos'] = A
population['Infectados'] = I
population['Hospitalizados'] = H
population['Graves'] = G
population['Recuperados'] = R
population['Vacunados'] = V
population['No Generan Inmunidad'] = FV
population['Fallecidos'] = D
population['Total'] = S+E+A+I+H+G+R+V+FV+D

which = st.sidebar.multiselect('Datos para visualizar:',['Susceptibles','Expuestos','Asintomáticos','Infectados',
'Hospitalizados','Graves','Recuperados','Vacunados','No Generan Inmunidad','Fallecidos','Total'],default=['Total',
'Recuperados','Fallecidos','Infectados','Expuestos'])
if st.selectbox('Tipo de Gráfico',('Linea','Area'))=='Area':
    st.area_chart(population[which])

else:
    st.line_chart(population[which])

    
#to json
# result = population.to_json(orient="records")
# parsed = json.loads(result)
# st.write(json.dumps(parsed, indent=4))
