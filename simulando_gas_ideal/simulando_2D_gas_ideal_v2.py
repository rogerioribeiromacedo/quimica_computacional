# -*- coding: utf-8 -*-
"""
Simulando gás ideal.

Algumas considerações:
    1) As partículas serão representadas como círculos de mesmo raio e mesma massa.
    2) As partículas serão posicionadas de maneira randômica e terão velocidades iniciais diferentes.
    3) O número de partículas dentro do espaço é fixo.
    4) As partículas sofre colisão entre si e com a paredes. As colisões são perfeitamente elásticas.

O programa segue algumas ideias:
    1) Existe um número fixo de partículas
    2) As partículas tem suas posições e velocidades iniciadas no t = 0
    3) A cada passo de interação:
        a) verifica se a partículo colide seja entre si ou com a parede
        b) atualiza a velocidade considerando conservação de energia e o momento linear
        b) nova posição da partícula é a posição anterior mais velocidade que é multiplicada pelo passo de tempo:
            - r(t + dt) = r(t) + v(t)dT

Sobre as colisões entre as partículas:
    1)

Sobre a colisão da partículas com a parede:
    1) O componente da velocidade que for perpendicular à parede é invertido.

Author.............: Rogério Ribeiro Macêdo.
Curriculum Lattes..: http://lattes.cnpq.br/8806221981552346
Última atualização.: 24 de Julho de 2024.
"""
# pylint: disable=invalid-name
# pylint: disable=import-error
try:
    import sys
    import numpy as np
    from itertools import product

    # Matplotlib
    import warnings
    warnings.filterwarnings("ignore")
    import matplotlib.pyplot as plt
    plt.rcParams['figure.dpi'] = 300
    from matplotlib.animation import FuncAnimation
except ImportError as e:
    print('[!] The required Python libraries could not be imported:', file=sys.stderr)
    print(f'\t{e}')
    sys.exit(1)

# Constantes
n_ljust = 50  # Tamanho do texto
k_b = 1.38*(10**(-23))


class GasIdeal:
    """Classe que descreve o gás ideal e a funções necessárias para extra."""

    def __init__(self, n_particulas, massa, raio, largura, v_inicial, duracao, n_passos):
        """Inicializa propriedades."""
        self.n_particulas = n_particulas  # número de partículas
        self.massa = massa  # massa da partícula (kg)
        self.raio = raio  # raio da partícula (m)
        self.largura = largura  # largura da caixa de simulação (m)
        self.duracao = duracao  # duração da simulação (s)
        self.n_passos = n_passos  # número de passos
        self.dt = duracao / n_passos  # timestep (s)
        self.v_inicial = v_inicial  # velocidade inicial da partícula (m/s)

        #
        # Posição das particulas
        #

        # Cria as posições das partículas em uma grade
        tam_grade = int(np.ceil(np.sqrt(n_particulas)))

        # Espaço entre as partículas
        espaco = largura / tam_grade

        # O valor x do centro da partícula
        x = np.linspace(raio + espaco / 2, largura - raio - espaco/2, tam_grade)

        # Gera um vetor com formato x, y
        pos = list(product(x, x))

        # Inicializa as posiçoes das particulas
        self.posicoes = np.array(pos[:n_particulas])

        theta = np.random.uniform(0, 2*np.pi, size=n_particulas)
        vx, vy = self.v_inicial * np.cos(theta), self.v_inicial * np.sin(theta)
        self.velocidades = np.stack((vx, vy), axis=1)

    def verifica_colisao(self):
        """
        Verifica se há colisão entre as partículas ou entre partícula e parede.

        Ao se verificar uma colisão, atualiza as velocidades considerando conservação
        de energia e momento linear.
        """
        # Nova posicao das partículas
        pos_nova = self.posicoes + self.velocidades * self.dt

        # Avalia as colições com as paredes
        # Eixo X
        self.velocidades[pos_nova[:, 0] < self.raio, 0] *= -1  # colisão contra a parede da esquerda (x)
        self.velocidades[pos_nova[:, 0] > self.largura-self.raio, 0] *= -1  # colisão com a parede da direita (x)
        # Exio Y
        self.velocidades[pos_nova[:, 1] < self.raio, 1] *= -1  # colisão com a parede inferior (y)
        self.velocidades[pos_nova[:, 1] > self.largura-self.raio, 1] *= -1  # colisão com a parede superior (y)

        # Avaliando a colisões entre as partículas
        for i in range(self.n_particulas):
            for j in range(i+1, self.n_particulas):
                # Calcula a norma ou módulo do vetor para saber se ouve colisão
                if np.linalg.norm(pos_nova[i] - pos_nova[j]) < 2*self.raio:
                    rdiff = self.posicoes[i] - self.posicoes[j]
                    vdiff = self.velocidades[i] - self.velocidades[j]

                    # atualiza velocidade da partícula i
                    self.velocidades[i] = self.velocidades[i] - rdiff.dot(vdiff)/rdiff.dot(rdiff)*rdiff

                    # atualiza velocidade da partícula j
                    self.velocidades[j] = self.velocidades[j] + rdiff.dot(vdiff)/rdiff.dot(rdiff)*rdiff

    def passo(self):
        """Calculando as posições."""
        self.verifica_colisao()
        self.posicoes += self.velocidades * self.dt

    def simular(self):
        """Simulando a movimentação de um gás ideal."""
        # Matriz de posições e velocidades (inicializando) para todos os passos
        #
        pos_simul = np.zeros((self.n_passos, self.n_particulas, 2))
        vel_simul = np.zeros((self.n_passos, self.n_particulas))

        for n in range(self.n_passos):
            pos_simul[n, :, :] = self.posicoes
            vel_simul[n, :] = np.linalg.norm(self.velocidades, axis=1)

            # Passo
            self.passo()

        return pos_simul, vel_simul

    def energia_cinetica_media(self, v: float) -> float:
        """
        Calcula a energia cinética média.

        Parameters
        ----------
        m : float
            Massa da partícula.
        v : []
            Lista de velocidades.

        Returns
        -------
        float
            Valor da energia cinética média.

        """
        media_quad_vel = np.sum(v**2) / self.n_particulas
        return 0.5*self.massa*media_quad_vel

    def MaxwellBoltzmann(self, v):
        """Distribuição de Maxwell."""
        # Energia cinética média
        KE_avg = 1/2*self.massa*np.sum(self.v_inicial**2)

        kT = KE_avg  # temperature
        sigma_sq = kT/self.massa  # spread of the distribution
        f = np.exp(-v**2/(2*sigma_sq)) * v/sigma_sq  # Maxwell-Boltzmann distribution
        return f


def head_msg():
    """
    Header message.

    Returns
    -------
    None.

    """
    print()
    print("-".center(79, "-"))
    print(f'{"|":<1} {"UNIFEI - Universidade Federal de Itajubá":^75} '
          f'{"|":>1}')
    print(f'{"|":<1} {"LaQC - Laboratório de Química Computacional":^75} '
          f'{"|":>1}')
    print(f'{"|":<1} {" ":^75} {"|":>1}')
    print(f'{"|":<1} {">>> use o comando [sair] para terminar o programa <<<":^75} '
          f'{"|":>1}')
    print(f'{"|":<1} {" ":^75} {"|":>1}')
    print("-".center(79, "-"))
    print(f'{"|":<1} {"Simulando partículas de um gás ideal.":^75} '
          f'{"|":>1}')
    print("-".center(79, "-"))
    print()


def tchau():
    """
    We are leaving and saying goodbye (tchau, inté!!).

    Returns
    -------
    None.

    """
    print("")
    print("-".center(79, "-"))
    print(f'{"|":<1} {"UNIFEI - Universidade Federal de Itajubá":^75} '
          f'{"|":>1}')
    print(f'{"|":<1} {"LaQC - Laboratório de Química Computacional":^75} '
          f'{"|":>1}')
    print(f'{"|":<1} {" ":^75} {"|":>1}')
    print(f'{"|":<1} {"Inté!!!!":^75} '
          f'{"|":>1}')
    print(f'{"|":<1} {" ":^75} {"|":>1}')
    print("-".center(79, "-"))
    print("")
    sys.exit()


def gerar_img_inicial(gas: GasIdeal):
    """
    Gera uma imagem inicial da distribuição das partículas.

    Parameters
    ----------
    gas : GasIdeal
        Classe de Gas Ideal.

    Returns
    -------
    None.

    """
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    #
    # Desenhando as partículas
    #
    for (x, y) in gas.posicoes:
        ax.add_artist(plt.Circle((x, y), gas.raio, color='r'))

    # Spines -> as linhas dos eixos
    for axis in ['bottom', 'left', 'top', 'right']:
        ax.spines[axis].set_linewidth(2.5)
        ax.spines[axis].set_color('0.2')

    # Configurando o eixo x
    ax.set_xlim(0, gas.largura)
    ax.xaxis.set_ticks_position('none')
    plt.xticks(visible=False)

    # Configurando o eixo y
    ax.set_ylim(0, gas.largura)
    ax.yaxis.set_ticks_position('none')
    plt.yticks(visible=False)

    ax.set_aspect('equal')

    # Salvando gráfico
    plt.savefig('img_inicial.png',
                bbox_inches='tight',
                dpi=300,
                facecolor=ax.get_facecolor())
    plt.cla()
    plt.close()


def animate_positions(frame, ax1, plt, positions, gas, speeds, v):
    """Função para animação."""
    ax1.clear()

    for i in range(gas.n_particulas):
        x, y = positions[frame, i, 0], positions[frame, i, 1]
        circle = plt.Circle((x, y), gas.raio, fill=True)
        ax1.add_artist(circle)

    ax1.set_xlabel('$x$', fontsize=15)
    ax1.set_ylabel('$y$', fontsize=15)
    ax1.set_title('Animação do Gás Ideal', fontsize=15)
    ax1.set_xlim(0, gas.largura)
    ax1.set_ylim(0, gas.largura)
    ax1.set_aspect('equal')
    ax1.set_xticks([])
    ax1.set_yticks([])

    plt.tight_layout()


def exibe_resultados(gas: GasIdeal, vel_simul: np.ndarray) -> None:
    """
    Exibindo resultados da simulação.

    Parameters
    ----------
    gas : GasIdeal
        Objeto da classe.
    vel_simul : numpy.ndarray
        Array contendo as velocidades das partículas.

    Returns
    -------
    None
        Nada.

    """
    # Verifica se a energia cinética total é conservada
    print(f' - Energia cinética total conversada: {np.sum(vel_simul[0]**2):.3f}, {np.sum(vel_simul[-1]**2):.3f}')

    Ek_media = gas.energia_cinetica_media(vel_simul)
    print(f" - Energia cinética média: {Ek_media:12.4E}")
    temp = (Ek_media/k_b)*(3/2)
    print(f" - Temperatura: {temp:8.2f} K")


def main():
    """
    Principal function.

    Returns
    -------
    None.

    """
    try:
        # Número de Partículas
        n_particulas = input("Número de partículas "
                             "[100]".ljust(n_ljust, ".") + ": ").strip()
        if n_particulas == "sair":
            tchau()
        else:
            if len(n_particulas) > 0:
                n_particulas = int(n_particulas)
            else:
                n_particulas = 100

        # Massa da partícula
        massa_particula = input("Massa da partícula (kg) "
                                "[5.31e-26 kg]".ljust(n_ljust, ".") + ": ").strip()
        if massa_particula == "sair":
            tchau()
        else:
            if len(massa_particula) > 0:
                massa_particula = float(massa_particula)
            else:
                massa_particula = 5.31*10**(-26)

        # Raio da Partícula
        raio_particula = input("Raio da partícula (m) "
                               "[0.3]".ljust(n_ljust, ".") + ": ").strip()
        if raio_particula == "sair":
            tchau()
        else:
            if len(raio_particula) > 0:
                raio_particula = float(raio_particula)
            else:
                raio_particula = 0.3

        # Tamanho da caixa de simulação
        l_caixa = input("Largura da caixa de simulação (m) "
                        "[20.0]".ljust(n_ljust, ".") + ": ").strip()
        if l_caixa == "sair":
            tchau()
        else:
            if len(l_caixa) > 0:
                l_caixa = float(l_caixa)
            else:
                l_caixa = 20.0

        # Velocidade inicial das partículas
        v_inicial = input("Velocidade inicial (m/s) "
                          "[2.0]".ljust(n_ljust, ".") + ": ").strip()
        if v_inicial == "sair":
            tchau()
        else:
            if len(v_inicial) > 0:
                v_inicial = float(v_inicial)
            else:
                v_inicial = 2.0

        # Duração da animação
        duracao = input("Duração da animação (s) "
                        "[10]".ljust(n_ljust, ".") + ": ").strip()
        if duracao == "sair":
            tchau()
        else:
            if len(duracao) > 0:
                duracao = int(duracao)
            else:
                duracao = 10

        # Número de passos
        n_passos = input("Número de passos (s) "
                         "[500]".ljust(n_ljust, ".") + ": ").strip()
        if n_passos == "sair":
            tchau()
        else:
            if len(n_passos) > 0:
                n_passos = int(n_passos)
            else:
                n_passos = 500

        # Inicializa a classe do gás
        print(" - Criando objeto do gás ideal.")
        gas = GasIdeal(n_particulas,
                       massa_particula,
                       raio_particula,
                       l_caixa,
                       v_inicial,
                       duracao,
                       n_passos)

        # Gerar imagem com as posições iniciais das partículas
        print(" - Salvando imagem com estrutura inicial.")
        gerar_img_inicial(gas)

        print(" - Simulando...")
        pos_simul, vel_simul = gas.simular()

        # Animando
        print(" - Gerando animação")
        v = np.linspace(0, 35, 500)
        fig, ax1 = plt.subplots(1, 1, figsize=(12, 6))

        # # Intervalo para animação
        interval = duracao*1e3 / n_passos
        animation = FuncAnimation(fig,
                                  animate_positions,
                                  frames=n_passos,
                                  interval=interval,
                                  fargs=[ax1, plt, pos_simul, gas, vel_simul, v])
        animation.save('gas_ideal.mp4', writer='ffmpeg', fps=30)

        #
        # Resultados
        #
        exibe_resultados(gas, vel_simul)

        print("")
        print(" - Fim da simulação.")

    except ValueError as erro:
        print(f" - Erro: {erro}")
        sys.exit(-1)


if __name__ == "__main__":
    # Show header message.
    head_msg()

    # Main
    main()
