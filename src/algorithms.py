import os
import Bio
import numpy as np
from tqdm import tqdm

from abc import abstractmethod
from collections import deque
from enum import Enum

def score(x, y):
    if x == "-" or y == "-":
        return -2
    elif x == y:
        return 1
    elif x != y:
        return -1

class AlignmentAlgorithm():

    class Node:
        """
        Class represents node required in reconstruct optimal alignment path
        """
        def __init__(self, firstSeq, secondSeq, i, j):
            self.firstSeq = firstSeq
            self.secondSeq = secondSeq
            self.i = i
            self.j = j
    
    class Arrow(Enum):
        """
        Enum for representing arrows
        """
        LEFT = 1
        DIAG = 2
        TOP = 3

    ARROW_MASK = (Arrow.LEFT, Arrow.DIAG, Arrow.TOP)
    GAP = "-"

    def __init__(self, firstSeq, secondSeq, scoreFunction = score):
        self.firstSeq = "-" + firstSeq
        self.secondSeq = "-" + secondSeq
        self.lenFirstSeq = len(self.firstSeq)
        self.lenSecondSeq = len(self.secondSeq)
        
        self.scoreMatrix = np.ndarray(shape=(self.lenSecondSeq, self.lenFirstSeq), dtype=np.int)
        self.arrowMatrix = np.ndarray(shape=(self.lenSecondSeq, self.lenFirstSeq), dtype=list)
        
        self.scoreFunction = scoreFunction

    @abstractmethod
    def _max_method(self, maxValue):
        pass
    
    @abstractmethod
    def _initialize_first_row_and_col(self):
        pass
    
    @abstractmethod
    def _ret_score(self):
        pass

    def _fillCell(self, i, j):
        """
        Function reposnsible for filling every cell based on dynamic programming rule
        """
        left = self.scoreMatrix[i, j - 1] + self.scoreFunction(self.GAP, self.GAP)
        diag = self.scoreMatrix[i - 1, j - 1] + self.scoreFunction(self.firstSeq[j], self.secondSeq[i])
        top = self.scoreMatrix[i - 1, j] + self.scoreFunction(self.GAP, self.GAP)
        
        directsArray = np.array([left, diag, top])
        maxScore = self._max_method(directsArray.max())
        directMask = directsArray == maxScore
        
        self.scoreMatrix[i, j] = maxScore
        self.arrowMatrix[i, j] = [arrow for mask, arrow in zip(directMask, self.ARROW_MASK) if mask == True]

    def _path_create(self, node, arrow):
        if arrow == self.Arrow.LEFT:
            firstSign = self.firstSeq[node.j]
            secondSign = self.GAP
            vector = (0, -1)
        elif arrow == self.Arrow.DIAG:
            firstSign = self.firstSeq[node.j]
            secondSign = self.secondSeq[node.i]
            vector = (-1, -1)
        elif arrow == self.Arrow.TOP:
            firstSign = self.GAP
            secondSign = self.secondSeq[node.i]
            vector = (-1, 0)
        nextNode = self.Node(firstSign + node.firstSeq, secondSign + node.secondSeq, node.i + vector[0], node.j + vector[1])
        return nextNode

    def calculate_score(self):
        """
        Calculate score matrix (with dynamic programming rule)
        """
        self.scoreMatrix[0, 0] = 0
        self._initialize_first_row_and_col()
        for i in tqdm(range(1, self.lenSecondSeq)):
            for j in range(1, self.lenFirstSeq):
                self._fillCell(i, j)
        return self._ret_score()

    @abstractmethod
    def best_alignments_(self, maxN):
        """
        Function prints <maxN> optimal alignments
        """
        pass

class WatermanSmithAlgorithm(AlignmentAlgorithm):
     
    def __init__(self, firstSeq, secondSeq, scoreFunction = score):
        super().__init__(firstSeq, secondSeq, scoreFunction)
    
    def _max_method(self, maxValue):
        return max(maxValue, 0)

    def _initialize_first_row_and_col(self):
        for i in range(1, self.lenSecondSeq):
            self.scoreMatrix[i, 0] = 0
        for j in range(1, self.lenFirstSeq):
            self.scoreMatrix[0, j] = 0

    def _ret_score(self):
        return self.scoreMatrix.max()

    def best_alignments_(self, maxN = 2):
        stack = deque()
        maxScore = self.scoreMatrix.max()
        rows, cols = np.where(self.scoreMatrix == maxScore)
        for row, col in zip(rows, cols):
            stack.append(self.Node("", "", row, col))
        i = 0
        while stack and i < maxN:
            actNode = stack.pop()
            if self.arrowMatrix[actNode.i, actNode.j]:
                for arrow in self.arrowMatrix[actNode.i, actNode.j]:
                    stack.append(self._path_create(actNode, arrow))
            if self.scoreMatrix[actNode.i, actNode.j] == 0:
                print(actNode.firstSeq)
                print(actNode.secondSeq)
                i = i + 1
                

class NeedlemanWunschAlgorithm(AlignmentAlgorithm):
    
    def __init__(self, firstSeq, secondSeq, scoreFunction = score):
        super().__init__(firstSeq, secondSeq, scoreFunction)
    
    def _max_method(self, maxValue):
        return maxValue
    
    def _initialize_first_row_and_col(self):
        for i in range(1, self.lenSecondSeq):
            self.scoreMatrix[i, 0] = i * self.scoreFunction(self.GAP, self.GAP)
        for j in range(1, self.lenFirstSeq):
            self.scoreMatrix[0, j] = j * self.scoreFunction(self.GAP, self.GAP)

    def _ret_score(self):
        return self.scoreMatrix[self.lenSecondSeq - 1, self.lenFirstSeq - 1]

    def best_alignments_(self, maxN = 2):
        stack = deque()
        stack.append(self.Node("", "", self.lenSecondSeq - 1, self.lenFirstSeq - 1))
        i = 0
        while stack and i < maxN:
            actNode = stack.pop()
            if self.arrowMatrix[actNode.i, actNode.j]:
                for arrow in self.arrowMatrix[actNode.i, actNode.j]:
                    stack.append(self._path_create(actNode, arrow))
            if actNode.i == 0 and actNode.j == 0:
                print(actNode.firstSeq)
                print(actNode.secondSeq)
                i = i + 1
  