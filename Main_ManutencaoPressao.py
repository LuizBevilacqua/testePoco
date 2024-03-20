"""
Main File to calculate __Manutenção de Pressão__ para poço cilindrico e reservatório circular
By: Endriw Rafael & Tiffany Franzoi
Avaliação de Formações
"""
import numpy as np
from scipy import special as sp
import Gavsteh_Func
import matplotlib.pyplot as plt


class PressureCalculator:
    def __init__(self, rD, reD, l, tD):
        self.rD = rD  # Raio do poço cilíndrico
        self.reD = reD  # Lista de distâncias reservatório-poço
        self.l = l
        self.tD = tD  # Array de tempo adimensional
        self.pD = []  # Lista para armazenar os perfis de pressão

    def calculate_pressure(self):
        for j in self.reD:
            # Função de pressão func(u) que será avaliada para cada tempo adimensional
            func = lambda u: (
                (sp.i0(j * np.sqrt(u)) * sp.k0(self.rD * np.sqrt(u)))
                - (sp.k0(j * np.sqrt(u)) * sp.i0(self.rD * np.sqrt(u)))
            ) / (
                (
                    (u ** (3 / 2))
                    * (
                        sp.i1(np.sqrt(u)) * sp.k0(j * np.sqrt(u))
                        + (sp.i0(j * np.sqrt(u)) * sp.k1(np.sqrt(u)))
                    )
                )
            )
            pD_reD = []  # Lista para armazenar os resultados da pressão

            for i in self.tD:
                # Chama a função gavsteh_param para calcular a pressão adimensional
                pD_reD.append(Gavsteh_Func.gavsteh_param(self.l, func, i))

            self.pD.append(
                pD_reD
            )  # Adiciona o perfil de pressão para esta distância reservatório-poço

    def plot_pressure_profiles(self):
        fig, ax = plt.subplots()
        for i in range(len(self.pD)):
            ax.plot(self.tD, self.pD[i], label=f"reD = {self.reD[i]}")
            ax.set_xscale("log")  # Gráfico log
            ax.set_ylim(3, 9)  # Define os limites do eixo y
            ax.grid(color="gray", linestyle="-.")
        plt.xlabel("tD")
        plt.ylabel("pD")
        plt.title(
            f"Gráfico de Pressão adimensional por Tempo adimensional, para reservatório circular com poço cilíndrico. rD={self.rD}"
        )
        plt.legend(framealpha=1)
        plt.yticks(
            np.arange(3, 9, step=0.25)
        )  # Define os intervalos das marcas no eixo y
        plt.show()


def main():
    rD = 1  # Raio do poço cilíndrico
    reD = [
        100,
        200,
        300,
        400,
        500,
        600,
        700,
        800,
        900,
        1000,
        1200,
        1400,
        1600,
        1800,
        2000,
        3000,
    ]  # Distâncias reservatório-poço
    l = 16
    tD = np.linspace(10e2, 10e6, 10000)  # Tempo adimensional

    # Cria uma instância da classe PressureCalculator
    pressure_calculator = PressureCalculator(rD, reD, l, tD)

    # Calcula os perfis de pressão
    pressure_calculator.calculate_pressure()

    # Plota os perfis de pressão
    pressure_calculator.plot_pressure_profiles()


if __name__ == "__main__":
    main()
