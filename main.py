## PART ONE ##
from Biosystem import *
from Part import *
from Rate import *
from Pulse import *
import matplotlib.pyplot as plt
import numpy as np


## PART TWO ##

# These two will be used to store the simulation results
T = None  # Time
Y = None  # Y axis

tMax = 30  # max simulation time
ktl = 4  # translation constant
ktx = 1  # transcription constant
beta_p = 0.2  # protein degradation constant
beta_m = 0.7  # mRna degradation constant

# Reporter variable (GFP)
beta_p2 = 0.4


## PART THREE ##

# Creating the biosystem
panamaUni = BioSystem()


## PART FOUR ##

# Defining variables AKA "compositors"
dmRNAdt = panamaUni.addCompositor('mRNA', 0)  # mRNA produced
dPdt = panamaUni.addCompositor('P', 0)  # protein produced
dGFPdt = panamaUni.addCompositor('GFP', 0)  # reporter

# Defining constants
panamaUni.addConstant('ktl', ktl)
panamaUni.addConstant('ktx', ktx)
panamaUni.addConstant('beta_p', beta_p)
panamaUni.addConstant('beta_m', beta_m)
panamaUni.addConstant('beta_p2', beta_p2)


## PART FIVE ##

# Defining the ordinary differential equations
ODE1 = Part('mRNA production', [dmRNAdt], [Rate('ktx - beta_m * mRNA')])
ODE2 = Part('Protein production', [dPdt], [Rate('ktl * mRNA - beta_p * P')])
ODE3 = Part('GFP production', [dGFPdt], [Rate('ktl * mRNA - beta_p2 * GFP')])

# Adding each of the parts
panamaUni.addPart(ODE1)
panamaUni.addPart(ODE2)
panamaUni.addPart(ODE3)


## PART SIX ##

(T, Y) = panamaUni.run([0, tMax])  # Running the simulation


## PART SEVEN ##

# Plotting graphs #
fig, axs = plt.subplots(2)
axs[0].set_title("Protein and reporter concentration")
axs[0].plot(T, Y[:, panamaUni.compositorIndex('P')], label="Proteína")
axs[0].plot(T, Y[:, panamaUni.compositorIndex('GFP')], label="Reportero")
axs[0].set(xlabel='Time (Min.)', ylabel='Concentration (fold)')
axs[0].legend(("Protein", "Reporter"))

# Correlation
axs[1].set_title("Methanol concentration correlation")
axs[1].plot(Y[:, panamaUni.compositorIndex('P')], Y[:, panamaUni.compositorIndex('GFP')],
            ":", label="Correlación Proteína vs. Reportero")
axs[1].set(xlabel='Protein', ylabel='Promoter')
fig.tight_layout()

plt.show()
