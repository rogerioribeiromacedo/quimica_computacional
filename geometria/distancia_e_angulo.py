# -*- coding: utf-8 -*-
"""
Avalia distâncias e ângulos dos átomos presentes em arquivo '.gro' (gromacs).

@author: Rogério Ribeiro Macêdo

"""
# pylint: disable=import-error
import sys
from pathlib import Path
import time
import numpy as np
import pandas as pd


TAM_TEXTO = 50
TAM_TEXTO_PROC = 35


def cabecalho():
    """
    Imprime cabeçalho do script.

    Returns
    -------
    None.

    """
    print()
    print("-".center(80, "-"))
    print(f'{"|":<1} {"Calcula distâncias e ângulos dos átomos":^76} '
          f'{"|":>1}')
    print(f'{"|":<1} {"de arquivos .gro do Gromacs":^76} '
          f'{"|":>1}')
    print(f'{"|":<1} {" ":^76} {"|":>1}')
    print(f'{"|":<1} {"obs: digite [sair] para encerrar o programa ":<76} '
          f'{"|":>1}')
    print(f'{"|":<1} {" ":^76} {"|":>1}')
    print("-".center(80, "-"))
    print()


def existe_arquivo(arquivo_gro):
    """
    Apenas verifica se o arquivo informado existe.

    Parameters
    ----------
    arquivo_gro : string
        String contendo o local onde encontra-se o arquivo.

    Returns
    -------
    Bool
        True, se arquivo existe; false, se arquivo não existe.

    """
    retorno = False
    path = Path(arquivo_gro)
    if path.is_file() and path.suffix == ".gro":
        retorno = True

    return retorno


def salvar_dataframe(df_dados, nome_arquivo, cols=[]):
    """
    Salva os dados do dataframe em um arquivo csv.

    Parameters
    ----------
    df_dados : Pandas dataframe
        Dataframe do Pandas contendo os dados a serem salvos.
    nome_arquivo : string
        Nome do arquivo.
    cols : array, opcional
        Colunas dos dados. Padrão é [].

    Returns
    -------
    status : bool
        True, nenhum erro ocorreu; False, erro ao salvar.

    """
    status = True
    nome_arquivo = nome_arquivo + ".csv"

    try:
        if len(cols) >= 1:
            df_dados.to_csv(nome_arquivo, columns=cols, index=False)
        else:
            df_dados.to_csv(nome_arquivo, index=False)
    except OSError as msg_erro:
        print(f"Erro ao salvar {nome_arquivo}! Erro: {msg_erro}")
        status = False
    except Exception as msg_erro:
        print(f"Erro ao salvar {nome_arquivo}! Erro: {msg_erro}")
        status = False

    return status


def criar_df_atomos(arquivo_gro):
    """
    Cria um dataframe (Pandas) contendo todos os átomo do arquivo.

    Parameters
    ----------
    arquivo_gro : string
        Nome/local do arquivo '.gro'.

    Returns
    -------
    df_atomos : dataframe Pandas
        Dataframe do Pandas contendo a lista de átomos do arquivo.

    """
    atomos = []
    colunas = ["numero_residuo", "nome_residuo", "nome_atomo",
               "numero_atomo", "x", "y", "z"]
    indicador = "/"

    # Listando átomos
    with open(arquivo_gro, "r") as f_arquivo:
        descricao = f_arquivo.readline().strip()
        print(" + Descrição".ljust(TAM_TEXTO_PROC, ".") + ": " + f"{descricao}")
        qtde_atomos = int(f_arquivo.readline().strip())
        print(" + Total de átomos".ljust(TAM_TEXTO_PROC, ".") + ": " + f"{qtde_atomos}")

        num_linhas = 0

        for linha in f_arquivo:
            # Sensação de tempo passando
            if indicador == "/":
                indicador = "\\"
            else:
                indicador = "/"
            print(" + Carregando átomos".ljust(TAM_TEXTO_PROC, ".") + f": {indicador}", end="\r")
            time.sleep(0.0005)

            num_linhas += 1
            if num_linhas <= qtde_atomos:
                dados_linha = []

                # número do resíduo (5 posições, integer)
                dados_linha.append(linha[0:5].strip())

                # nome do resíduo (5 caracteres)
                dados_linha.append(linha[5:10].strip())

                # atom name (5 characters)
                dados_linha.append(linha[10:15].strip())

                # atom number (5 positions, integer)
                dados_linha.append(linha[15:20].strip())

                # position (in nm, x y z in 3 columns, each 8 positions
                # with 3 decimal places)
                dados_linha.append(float(linha[20:28].strip()))
                dados_linha.append(float(linha[28:36].strip()))
                dados_linha.append(float(linha[36:44].strip()))

                # velocity (in nm/ps (or km/s), x y z in 3 columns,
                # each 8
                # positions with 4 decimal places)
                # considering that the file don't have it

                atomos.append(dados_linha)
        f_arquivo.close()

    # Criando o dataframe
    df_atomos = pd.DataFrame(
        atomos,
        columns=colunas)

    # Salva dados
    print("")
    salvar_dataframe(df_atomos, "lista_atomos", colunas)
    print(" + Lista de átomos salva!")

    return df_atomos


def dist_oxi_oxi(df_from_gro):
    """
    Calcula distância entre os átomos de oxigênio.

    Parameters
    ----------
    df_from_gro : Pandas dataframe
        Dataframe contendo todos os átomos do sistema.

    Returns
    -------
    df_dist_oxi_oxi : Pandas dataframe
        Dataframe contendo as distâncias entre os átomos de oxigênio.

    """
    # colunas
    cols = ["oxigenio_A", "oxigenio_B", "distancia"]

    # filtro (lista todos os átomos de oxigênio)
    df_oxigens = df_from_gro[df_from_gro['nome_atomo'] == "OW"]

    # data
    data = []
    indicador = "/"

    # calculando distâncias
    for ind1 in df_oxigens.index:
        for ind2 in df_oxigens.index[ind1+1:]:
            # Sensação de tempo passando
            if indicador == "/":
                indicador = "\\"
            else:
                indicador = "/"
            print(" + Distância entre oxigênios".ljust(TAM_TEXTO_PROC, ".") + f": {indicador}", end="\r")
            time.sleep(0.0005)

            vlr_x2 = (df_oxigens['x'][ind1] - df_oxigens['x'][ind2]) ** 2
            vlr_y2 = (df_oxigens['y'][ind1] - df_oxigens['y'][ind2]) ** 2
            vlr_z2 = (df_oxigens['z'][ind1] - df_oxigens['z'][ind2]) ** 2

            distance = np.sqrt((vlr_x2 + vlr_y2 + vlr_z2))

            data.append([df_oxigens["nome_atomo"][ind1] + "_" +
                         df_oxigens["numero_residuo"][ind1],
                         df_oxigens["nome_atomo"][ind2] + "_" +
                         df_oxigens["numero_residuo"][ind2],
                         np.round(distance, 2)])

    # Criando dataframe
    df_dist_oxi_oxi = pd.DataFrame(
        data,
        columns=cols)

    # Salva dados
    print("")
    salvar_dataframe(df_dist_oxi_oxi, "dist_oxi_oxi", cols)
    print(" + Distâncias entre oxigênios salva!")

    return df_dist_oxi_oxi


def calc_angle(ponto1, ponto2, ponto3):
    """
    Calcula o ângulo entre três pontos.

    Parameters
    ----------
    ponto1 : TYPE
        DESCRIPTION.
    ponto2 : TYPE
        DESCRIPTION.
    ponto3 : TYPE
        DESCRIPTION.

    Returns
    -------
    vetor1 : TYPE
        DESCRIPTION.
    vetor2 : TYPE
        DESCRIPTION.
    angle : TYPE
        DESCRIPTION.

    """
    ponto1 = np.array(ponto1)
    ponto2 = np.array(ponto2)
    ponto3 = np.array(ponto3)

    vetor1 = ponto2 - ponto1
    vetor1 = np.around(vetor1, 4)
    vetor2 = ponto2 - ponto1
    vetor2 = np.around(vetor2, 4)

    norma_vetor1 = np.sqrt(vetor1[0]**2 + vetor1[1]**2 + vetor1[2]**2)
    norma_vetor1 = np.around(norma_vetor1, 4)
    norma_vetor2 = np.sqrt(vetor2[0]**2 + vetor2[1]**2 + vetor2[2]**2)
    norma_vetor1 = np.around(norma_vetor1, 4)

    produto_vetorial = vetor1.dot(vetor2)

    cos_angle_rad = produto_vetorial / (norma_vetor1 + norma_vetor2)
    angle = np.degrees(np.arccos(cos_angle_rad))
    angle = np.around(angle, 4)

    return vetor1, vetor2, angle


def molecules_angles(df_from_gro):
    """
    Calcula angulo entre as moléculas de água.

    Parameters
    ----------
    df_from_gro : TYPE
        DESCRIPTION.

    Returns
    -------
    df_molecules_angles : TYPE
        DESCRIPTION.

    """
    cols = ["numero_residuo", "vetor_1", "vetor_2", "angulo"]

    # filtro (pega todos os átomos de oxigênio [OW])
    df_all_ow = df_from_gro[df_from_gro['nome_atomo'] == "OW"]

    # filtro (pega todos os átomos de hidrogêncio [HW] )
    df_all_hw = df_from_gro[df_from_gro['nome_atomo'].str.startswith("HW")]

    # data
    data = []
    indicador = "/"

    for ind1 in df_all_ow.index:
        # Sensação de tempo passando
        if indicador == "/":
            indicador = "\\"
        else:
            indicador = "/"
        print(" + Ângulos das moléculas de água".ljust(TAM_TEXTO_PROC, ".") + f": {indicador}", end="\r")
        time.sleep(0.0005)

        residue_number = df_all_ow["numero_residuo"][ind1]
        ponto_origem = [df_all_ow["x"][ind1], df_all_ow["y"][ind1], df_all_ow["z"][ind1]]

        pontos_destino = []
        for ind2 in df_all_hw[df_all_hw["numero_residuo"] == residue_number].index:
            pontos_destino.append([df_all_hw["x"][ind2], df_all_hw["y"][ind2], df_all_hw["z"][ind2]])

        angle = calc_angle(ponto_origem, pontos_destino[0], pontos_destino[1])
        data.append([residue_number, angle[0], angle[1], angle[2]])

    # Criando dataframe
    df_molecules_angles = pd.DataFrame(
        data,
        columns=cols)

    # Salva dados
    print("")
    salvar_dataframe(df_molecules_angles, "molecules_angles", cols)
    print(" + Ângulos das moléculas salva!")

    return df_molecules_angles


def main(arquivo_gro):
    """
    Procedimento principal.

    Parameters
    ----------
    arquivo_gro : string
        Local do arquivo .gro.

    Returns
    -------
    None.

    """
    df_atomos = criar_df_atomos(arquivo_gro)

    # Calculando a distância entre os átomos de oxigênio
    df_dist_oxi_oxi = dist_oxi_oxi(df_atomos)
    oxi_oxi_mean = df_dist_oxi_oxi["distancia"].mean()

    # calculando ângulos (água)
    df_molecules_angles = molecules_angles(df_atomos)
    angle_mean = df_molecules_angles["angulo"].mean()
    angle_max = df_molecules_angles["angulo"].max()
    angle_min = df_molecules_angles["angulo"].min()

    # Imprimindo o resumo
    print("")
    print("-".center(80, "-"))
    print(f'{"|":<1} {"Resumo":<76} '
          f'{"|":>1}')
    print("-".center(80, "-"))
    print("Média de distância entre os oxigênios".ljust(40, ".") + ": " + f" {oxi_oxi_mean:1.6} angstrons")
    print("Mean value of angles".ljust(40, ".") + ": " + f" {angle_mean:1.6} degrees")
    print("Max value of angles".ljust(40, ".") + ": " + f" {angle_max:1.6} degrees")
    print("Min value of angles".ljust(40, ".") + ": " + f" {angle_min:1.6} degrees")


if __name__ == '__main__':
    cabecalho()

    if len(sys.argv) == 2:
        if existe_arquivo(sys.argv[1]):
            main(sys.argv[1])
        else:
            print(f" + Arquivo ({sys.argv[1]}) não existe!")
            sys.exit()
    else:
        arquivo = input("Local e nome do arquivo de estrutura "
                        "(.gro)".ljust(TAM_TEXTO, ".") + ": ").strip()
        if existe_arquivo(arquivo):
            main(arquivo)
        else:
            print(f" + Arquivo ({arquivo}) não existe!")
            sys.exit()
