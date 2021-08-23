from collections import deque
import numpy as np
import ntpath
import os
from time import time

def combine_columns(Spart):
    '''
    Funzione che combina due o più colonne nel calcolo del vettore rappresentativo di un insieme lamda
    composto da due o più elementi del dominio
    :param singletonMatrixPartition: partizione della matrice dei vettori rappresentativi dei vettori singoletti
    :return: il vettore rappresentativo dell'insieme lamda composto da due o più elementi di M
    '''
    res = [0]*Spart.shape[0]
    while Spart.shape[1] >= 2:
        c1 = Spart[:, 0]
        c2 = Spart[:, 1]
        for i in range(Spart.shape[0]):
            if c1[i] != 0 and c2[i] != 0:
                res[i] = -1
            elif c1[i] == -1 or c2[i] == -1:
                res[i] = -1
            elif c1[i] != 0 and c2[i] == 0:
                res[i] = c1[i]
            elif c1[i] == 0 and c2[i] != 0:
                res[i] = c2[i]
            #elif c1[i] == 0 and c2[i] == 0:
                #res[i] = 0
        Spart = Spart[:, 2:]
        Spart = np.concatenate((np.array(res).reshape(len(res),1), Spart), axis=1)
    return Spart.flatten().tolist()

def build_projection(lamda, representativeVector):
    '''
    Costruisce la projection del contenuto del
    vettore rappresentativo sull'insieme lamda associato
    :param lamda: insieme associato al vettore rappresentativo
    :param representativeVector: vettore rappresentativo di lamda
    :return: la projection intesa come insieme di elementi di lamda
    contenuti nel vettore rappresentativo
    '''
    projection = []
    for elem in representativeVector:
        if elem in lamda:
            projection.append(elem)
    return set(projection)
    

def build_representativeVector(lamda, S):
    '''
    Costruisce il vettore rappresentativo di un insieme lamda
    :param lamda: insieme di cui costruire il vettore rappresentative
    :param singletonRepresentativeMatrix: matrice dei vettori rappresentativi dei
    sottoinsiemi singoletti di M
    :return: il vettore rappresentativo
    '''
    if len(lamda)==0:
        return [0]*S.shape[0]
    elif len(lamda)==1:
        return np.array(S, dtype=np.int64)[:, lamda[0]-1].tolist()
    else:
        return combine_columns(np.array(S, dtype=np.int64)[:, [x-1 for x in lamda]])

def check(lamda, S):
    '''
    Effettua il controllo  CHECK sull'insieme lamda che gli viene passato
    :param lamda: insieme lamda di cui effettuare il CHECK
    :param singletonRepresentativeMatrix matrice dei vettori rappresentative dei singoletti
    :return: OK se l'insieme è un potenziale hitting set, KO se l'insieme non è un
    hitting set né può diventarlo, MHS se l'insieme è già un minimum hitting set
    '''
    representativeVector = build_representativeVector(lamda, S)
    projection = build_projection(lamda, representativeVector)

    if projection == set(lamda):
        if np.count_nonzero(representativeVector) == len(representativeVector):
            return 'MHS'
        else:
            return 'OK'
    else:
        return 'KO'
    
def output(lamda, countMHS, mapping):
    '''
    effettua l'ouput dell'insieme lamda, rivelatosi un mhs, insieme alla sua cardinalità e al conteggio corrente di mhs trovati
    :param lamda: di cui effettuare l'output
    :return: l'output di lamda mhs
    '''
    if mapping is None:
        print('MHS found: {} of dimension {}'.format(lamda, len(lamda)))
        print('MHS encountered : {} \n'.format(countMHS))
    else:
        lamda_remapped = [mapping[elem-1] for elem in lamda]
        print('MHS found: {} of dimension {}'.format(lamda_remapped, len(lamda_remapped)))
        print('MHS encountered : {} \n'.format(countMHS))

def getMaps(indecesRemoved, MprimeLength):
    '''
    Funzione che calcola una mappa delle colonne di A dal dominio M al dominio M'
    quando si è effettuata una pre-elaborazione
    :param indecesRemoved: indici delle colonne rimosse nelle pre-elaborazione
    :param MprimeLength: lunghezza del nuovo dominio M'
    :return:
    '''
    indecesRemoved = [idx + 1 for idx in indecesRemoved]
    acc = 0
    mapping = []
    i = 1
    while len(mapping) < MprimeLength:
        while acc + i in indecesRemoved:
            acc += 1
        mapping.append(i+acc)
        i += 1
    return mapping

def getSingletonRepresentativeMatrix(A):
    '''
    Restituisce una matrice, intesa come vettore di vettori, contenente i vettori
    rappresentativi degli sotto-insiemi singoletti di M
    :param A: matrice di input
    :return: i vettori rappresentativi degli insiemi singoletti sottoforma di matrice
    '''
    for i in range(len(A)):
        for j in range(len(A[0])):
            if A[i][j] == 1:
                A[i][j] = j + 1
    return A

def del_rows(A):
    toBeRemoved = []
    i = 0
    while i <= A.shape[0]-2:
        j = i+1
        while j <= A.shape[0]-1:
            diff = A[i] - A[j]
            diff2 = A[j] - A[i]
            if np.count_nonzero(diff == -1) == 0:
                toBeRemoved.append(i)
            elif np.count_nonzero(diff2 == -1) == 0:
                toBeRemoved.append(j)
            j += 1
        i += 1
    print('Rows dropped in pre-processing: {} \n'.format([x+1 for x in toBeRemoved]))
    return np.delete(A, toBeRemoved, axis=0)

def del_cols(A):
    zeroCols =(~np.all(A==0, axis=0)).tolist()
    indecesRemoved = []
    for i in range(len(zeroCols)):
        if zeroCols[i] == False:
            indecesRemoved.append(i)
    print('Columns dropped in preprocessing: {} \n'.format([idx+1 for idx in indecesRemoved]))
    return (A[:, zeroCols], indecesRemoved)

def pre_processing(A):
    '''
    Restituisce una tupla contentene la matrice A di N' righe e M' colonne a seguito della
    pre-elaborazione e gli indici delle colonne che sono state rimosse
    :param A: matrice A input dell'algoritmo
    '''
    start = time()
    (Aprime, indecesRemoved) = del_cols(del_rows(np.array(A, dtype=np.int64)))
    end = time()
    print("Preprocessing required %.6f seconds to execute \n" % (end - start))
    return (Aprime, indecesRemoved)

def mbase(A, timeEnabled=True, mapping=None):
    '''
    Algoritmo mbase per il calcolo di tutti i mhs data una matrice A che riassume il contenuto del dominio M
    e delle collezione di insieme N
    :param A: matrice di ingresso sopracitata
    :param timeEnabled: flag per il calcolo del tempo necessario alla computazione
    :param mapping: mapping delle colonne di A qualora si sia compiuta una pre-elaborazione
    '''
    if timeEnabled:
        start = time()

    coda = deque()
    coda.append([])
    countMHS=0

    S = getSingletonRepresentativeMatrix(np.copy(A))
    M = list(range(1, len(A[0])+1))

    while len(coda)>0:
        alpha = coda.popleft()

        if len(alpha)==0:
            e = M[0]
        else:
            e = max(alpha)+1

        while e <= M[-1]:
            lamda = alpha + [e]
            result = check(lamda, S)

            if result == 'OK' and e != M[-1]:
                coda.append(lamda)
            elif result == 'MHS':
                countMHS += 1
                output(lamda, countMHS, mapping)
                
            e += 1
    if timeEnabled:
        end = time()
        print("MBASE required %.6f seconds to execute \n Execution completed \n" % (end-start))


def getMatrixFromFile(filename):
    '''
    Legge un file .matrix dato il filename e restituisce una matrice a partire
    dal suo contenuto
    Per esempio a partire da './benchmarks1/74L85.000.matrix' si analizza il
    contenuto del file 74L85.000.matrix per creare la matrice
    :param filename: nome del file .matrix (indirizzo completo)
    :return: la matrice A costruita a partire dal contenuto del file
    '''
    print(ntpath.basename(filename))
    if os.stat(filename).st_size == 0:
        return np.array([[]],dtype=np.int64)
    else:
        file = open(filename, 'r')
        lines = file.readlines()
        lines = lines[5:]
        lines = [str.replace(line, ' -\n', '') for line in lines]
        lines = [str.replace(line, ' ', ',') for line in lines]
        #adesso abbiamo le righe della matrice A e possiamo costruirla
        rows = []
        for line in lines:
            rows.append(np.fromstring(line, sep=','))
        A = np.array(rows, dtype=np.int64)
        file.close()
        return A.tolist()

A = getMatrixFromFile(filename=input("Insert filename: \n"))
if len(A) == 0:
    print("The specified file is empty. Can't start the computation")
else:
    print('(|N| = {}, |M| = {}) \n'.format(len(A), len(A[0])))
    mbase(A)
    
    (A_post, indecesRemoved) = pre_processing(A)
    print('(|N| = {}, |M| = {}) \n'.format(len(A_post), len(A_post[0])))
    mbase(A_post, mapping=getMaps(indecesRemoved, len(A_post[0])))
