import time
import random

from typing import Iterable
from settings import *
from settings import Settings
from unionfind import DisjointSet


import graph

CYKLES_LIMIT_REACHED: int = 100
PATH_FOUND          : int = 101
SIMULATION_RUNNING  : int = 102

class Symulacja2D:
    NS = 0
    WE = 1

    def __init__(self, settings: Settings) -> None:
        '''
        Fields:
        ----------
        __frame: list[list[Cell]]
            list (matrix) of Cells

        __changed: list[tuple[int, int, int]]
            list od cells that have changed since previoues generation

        __width: int
            world width ( number of cells in axis X)

        __height: int
            world height (number of cells in axix Y)

        __datacollector: Datacollector
            for data collection

        __unionfind: UnionFind
            for finding when defects accumulated and destroyed grid


        __found_path: list(tuples[int, int, int])
            contains list of cells that destroyed grid
     
        '''
        self.__frame: list[list[int]] = []
        self.__changed: list[tuple(int, int, int)] = []
        self.settings: Settings = settings

        self._width, self._height  = self.settings.GRID_SIZE

        self._disjoint_set_AC = DisjointSet()
        self._disjoint_set_BD = DisjointSet()
        
        self.edge_cells = set()
        self._seed = self._width
        self._path_starting_point: int = None

        self.make_world()
        self._make_disjoint_sets()

        self._P1_counter: int = self._width * self._height
        self._P2_counter: int = 0
        self._PD_counter: int = 0
        self._D_counter : int = 0

        self.simulation_running: bool = True
        self.cykles_counter: int = 0

        self.connectivity_test_result = None
        self.defecty_d = []
        
    
    def __str__(self):
        frame1: list[list[int]] = self.__frame
        s1 = '\n'.join([''.join(['{:6}'.format(item) for item in row]) for row in frame1])

        return '\n' + s1 + '\n'


    def _encode_id(self, y, x):
        return y* (self._seed + 2) + x


    def _decode_id(self, id_: int):
        x, y = divmod(id_, self._seed + 2)
        return x, y
        

    @staticmethod
    def _vicinity(y:int, x:int) -> tuple[int, int]:
        '''
        Return coordinates of neighbor cell (von Neumann definition)
        
        Parameters:
        ----------
        x: int
            x - axis coordinate of base cell

        y: int
            y - axis coordinate of base cell

        Returns:
        -------
            tuple(int, int) - calculated coordinates of neighbores cells
        '''

        yield (y + 1, x) # Top neighbor
        yield (y - 1, x) # Bottom neighbor
        yield (y, x + 1) # Right neighbor
        yield (y, x - 1) # Left neighbor


  

    def _make_disjoint_sets(self) -> None:
        self.__a = set()
        self.__b = set()
        self.__c = set()
        self.__d = set()
        encode = self._encode_id

        w = self._width
        h = self._height

        unionAC = self._disjoint_set_AC
        unionBD = self._disjoint_set_BD

        unionAC.makeset(self.settings.SUPERCELL_A)
        unionBD.makeset(self.settings.SUPERCELL_B)
        unionAC.makeset(self.settings.SUPERCELL_C)
        unionBD.makeset(self.settings.SUPERCELL_D)

        for i in range(1, w + 1):
            #A
            cell_id = encode(0, i)
            unionAC.makeset(cell_id)
            unionAC.union(cell_id, self.settings.SUPERCELL_A)
            self.__a.add(cell_id)

            #C
            cell_id = encode(h + 1, i) 
            unionAC.makeset(cell_id)
            unionAC.union(cell_id, self.settings.SUPERCELL_C)
            self.__c.add(cell_id)

        for i in range(1, h + 1):
            #B
            cell_id = encode(i, w + 1)
            unionBD.makeset(cell_id)
            unionBD.union(cell_id, self.settings.SUPERCELL_B)
            self.__b.add(cell_id)


            #D
            cell_id = encode(i, 0)
            unionBD.makeset(cell_id)
            unionBD.union(cell_id,self.settings.SUPERCELL_D)
            self.__d.add(cell_id)
            
        self.edge_cells  = self.__a | self.__b | self.__c | self.__d



    @property
    def cykles(self) -> int:
        return self.cykles_counter


    @property
    def PD(self):
        return self._PD_counter

    @property
    def P1(self):
        return self._P1_counter


    @property
    def P2(self):
        return self._P2_counter

    
    @property
    def D(self):
        return self._D_counter

    @property
    def paths_direction(self):
        return self.connectivity_test_result



    def simulation_finished(self) -> bool:
        '''
        Checks if simulation ended. Simulation ends when defects riches oposite edges

        Returns:
        -------
        bool
            True if defects damaged grid otherwise False
        '''
        A = self.settings.SUPERCELL_A
        B = self.settings.SUPERCELL_B
        C = self.settings.SUPERCELL_C
        D = self.settings.SUPERCELL_D

        AC = self._disjoint_set_AC
        BD = self._disjoint_set_BD

        testNS = AC.connected(A, C)
        testWE = BD.connected(B, D)

        if testNS or testWE:
            self.connectivity_test_result = (testNS, testWE)
            if testNS and testWE:
                choice = random.choice([Symulacja2D.NS, Symulacja2D.WE])
                if choice == Symulacja2D.NS:
                    self._path_starting_point = A
                else:
                    self._path_starting_point = B

            elif testWE:
                self._path_starting_point = B

            elif testNS:
                self._path_starting_point = A

            return True

        else:
             return False

       

    def simulation_state(self) -> int:
        # Sprawdź czy znaleziono ścieżke perkolacyjną
        if self.simulation_finished():
            return PATH_FOUND

        # wtakim razie symulacja trwa dalej
        return SIMULATION_RUNNING



    def get_fracture_path(self, mark_connected_nodes: bool = False) -> Iterable[tuple[int, int, int]]:

        start = time.time()
        A: int = self.settings.SUPERCELL_A
        B: int = self.settings.SUPERCELL_B
        C: int = self.settings.SUPERCELL_C
        D: int = self.settings.SUPERCELL_D

        AC: DisjointSet = self._disjoint_set_AC
        BD: DisjointSet = self._disjoint_set_BD

        edge_A: set = self.__a
        edge_B: set = self.__b
        edge_C: set = self.__c
        edge_D: set = self.__d

        supercells = set([A,B,C,D])

        Graph = graph.Graph(weighted=False, digraph=False)
        starting_point: int = self._path_starting_point
        cells = []
        path  = []


        #wez wszystkie weżły które utworzyły rozdarcie w siatce
        if starting_point == A or starting_point == C:
            cells = AC.return_members_set(starting_point)

        elif starting_point == B or starting_point == D:
            cells = BD.return_members_set(starting_point)

        # print("PATHFINDING: Building graph -> start")
        building_graph = time.time()
        for cell in cells:
            Graph.add_node(cell)


        for node in Graph.nodes():
            y, x = self._decode_id(node)
            for b, a in self._vicinity(y, x):

                neighbor = self._encode_id(b, a)
                if neighbor in Graph \
                    and not (neighbor in edge_B and node in edge_D) \
                    and not (neighbor in edge_D and node in edge_B):
                    Graph.add_edge(node, neighbor)

        #połącz superkomórki z komórkami brzegowymi
        if starting_point == A or starting_point == C:
            for node in edge_A:
                if node in Graph:
                    Graph.add_edge(A, node)

            for node in edge_C:
                if node in Graph:
                    Graph.add_edge(C, node)

            path = Graph.shortest_path_bfs(C, A)

        elif starting_point == B or starting_point == D:
            for node in edge_B:
                if node in Graph:
                    Graph.add_edge(B, node)

            for node in edge_D:
                if node in Graph:
                    Graph.add_edge(D, node)

            path = Graph.shortest_path_bfs(D, B)
        
        bg = f'PATHFINDING: Building graph -> end {time.time() - building_graph: 6f}'
        # print(bg)

        returned_elements = []

        #usun superkomorki i komorki brzegowe z listy
        if mark_connected_nodes:
            for elem in cells:
                if elem not in supercells and elem not in self.edge_cells:
                    y, x = self._decode_id(elem)
                    returned_elements.append((y - 1, x - 1, self.settings.DESTRUCTION_PATH))

        for elem in path:
            if elem not in supercells and elem not in self.edge_cells:
                y, x = self._decode_id(elem)
                returned_elements.append((y - 1, x - 1, self.settings.FRACTURE_PATH))

        s = f'PATHFINDING: All {time.time() - start: 6f}'
        # print(s)

        return returned_elements

    

    def make_world(self):
        # print("SIMULATION: BUILDING GRID")
        start = time.time()
        yaxis = self._height + 2
        xaxis = self._width + 2
        frame = self.__frame

        for row in range(yaxis):
            line = []
            for col in range(xaxis):
                if row > 0 and col > 0 and col < self._width + 1 and row < self._height + 1:    
                    line.append(self.settings.P1)

                else:
                    line.append(self.settings.EDGE)

            frame.append(line)

        s = f'\n SUMULATION: GRID BUILDING -> end {time.time() - start}'
        # print(s)
                
               

    def _changed(self, y, x, state) -> None:
        '''
        co

        Parameters:
        ----------
        x: Cell
            Cell object from which we wll extract infomrmation about location
            and category of a cell
        '''

        self.__changed.append(( y - 1, x - 1, state))

    @staticmethod
    def neighbors(y, x):
        yield (y - 1, x)
        yield (y + 1, x)
        yield (y, x - 1)
        yield (y, x + 1)


    def neighbors2(self, y, x, line, next_frame):
        yield (y - 1, x, next_frame[y - 1][x])
        yield (y + 1, x, self.__frame[y + 1][x])
        yield (y, x - 1, line[x - 1])
        yield (y, x + 1, self.__frame[y][x + 1])



    def nerby_defect(self, y, x):
        f = self.__frame
        s = 0
        DEFECT = self.settings.DEFECT

        for b, a in self.neighbors(y, x):
            if f[b][a] < 0:
                s = s + DEFECT
            elif f[b][a] == DEFECT:
                s = s + 0
            else:
                s = s + f[b][a]

        return True if s > DEFECT else False




    def losuj(self, marker_P1, marker_P2):
        frame = self.__frame

        s = self.settings

        for y in range(1, self._height + 1):
            for x in range( 1, self._width + 1):
                cell = self.__frame[y][x]

                if cell == s.P1:
                    if random.random() < s.P1_PROBABILITY:
                        frame[y][x] = marker_P1

                
                elif cell == s.P2:
                    if random.random() < s.P2_PROBABILITY:
                        frame[y][x] = marker_P2

                elif cell > s.DEFECT and cell < s.PERMANENT:
                    frame[y][x] = frame[y][x] - 1
                

  



    def new_cell_state(self, markerp1, markerp2):

        s = self.settings
        next_frame = []

        for y in range(self._height + 2):
            line = []
            for x in range(self._width + 2):
                cell = self.__frame[y][x]
                if x > 0 and y > 0 and x < self._width + 1 and y < self._height + 1:
                    
                    if cell == markerp1 or cell == markerp2:
                        if self.nerby_defect(y, x) == True:
                            cell = s.PERMANENT
                            self._changed(y, x, s.PERMANENT)
                            enc = self._encode_id(y, x)
                            self._disjoint_set_AC.makeset(enc)
                            self._disjoint_set_BD.makeset(enc)

                            if line[x - 1] == s.P1:
                                line[x -1] = s.P2
                                self._changed(y, x - 1, s.P2)
                            
                            if next_frame[y - 1][x] == s.P1:
                                next_frame[y - 1][x] = s.P2
                                self._changed(y - 1, x, s.P2)

                            for b, a, state in self.neighbors2(y, x, line, next_frame):
                                n = state
                                if n == s.PERMANENT or n == s.EDGE:
                                    nid = self._encode_id(b, a)
                                    self._disjoint_set_AC.union(enc, nid)
                                    self._disjoint_set_BD.union(enc, nid)

                            if cell == markerp1:
                                self._P1_counter -= 1
                            else:
                                self._P2_counter -= 1
                            self._PD_counter += 1

                        else:
                            cell = s.DEFECT + s.HEAL_CYKLES
                            if line[x - 1] == s.P1:
                                line[x - 1] = s.P2
                                self._changed(y, x - 1, s.P2)

                            if next_frame[y - 1][x] == s.P1:
                                next_frame[y - 1][x] = s.P2
                                self._changed(y - 1, x, s.P2)
                            
                            # if self.__frame[y][x + 1] == s.P1:
                            #     self.__frame[y][x+1] = s.P2
                            #     self._changed(y, x + 1, s.P2)

                            # if self.__frame[y+1][x] == s.P1:
                            #     self.__frame[y+1][x] = s.P2
                            #     self._changed(y + 1, x, s.P2)
                            
                            self._changed(y, x, s.DEFECT)
                            if cell == markerp1:
                                self._P1_counter -= 1
                            else:
                                self._P2_counter -= 1
                            self._D_counter += 1


                    elif cell == s.P1:
                        if self.nerby_defect(y, x) == True:
                            cell = s.P2
                            self._changed(y, x, cell)
                            self._P1_counter -= 1
                            self._P2_counter += 1

                        else:
                            cell == s.P1

                    elif cell == s.P2:
                        if self.nerby_defect(y, x) == True:
                            cell = s.P2

                        else:
                            cell = s.P1
                            self._changed(y, x, cell)
                            self._P2_counter -= 1
                            self._P1_counter += 1

                    elif cell == s.DEFECT:
                        if line[x - 1] >= s.DEFECT or next_frame[y - 1][x] >= s.DEFECT or \
                            self.__frame[y + 1][x] >= s.DEFECT or self.__frame[y][x + 1] >= s.DEFECT:
                            cell = s.P2
                        else:
                            cell = s.P1
                        self._changed(y, x, cell)
                        self._P1_counter += 1
                        self._D_counter -= 1
                      

                    elif cell > s.DEFECT and cell < s.PERMANENT:
                        if self.nerby_defect(y, x) == True:
                            cell = s.PERMANENT
                            self._changed(y, x, cell)
                            # if line[x - 1] == s.P1:
                            #     line[x - 1] = s.P2
                            #     self._changed(y, x - 1, s.P2)

                            # if next_frame[y-1][x]==s.P1:
                            #     next_frame[y-1][x] = s.P2
                            #     self._changed(y-1, x, s.P2)

                            enc = self._encode_id(y, x)
                            self._disjoint_set_AC.makeset(enc)
                            self._disjoint_set_BD.makeset(enc)

                            for b, a, state in self.neighbors2(y, x, line, next_frame):
                                n = state
                                if n == s.PERMANENT or n == s.EDGE:
                                    nid = self._encode_id(b, a)
                                    self._disjoint_set_AC.union(enc, nid)
                                    self._disjoint_set_BD.union(enc, nid)
                            
                            self._D_counter -= 1
                            self._PD_counter += 1

                        

                line.append(cell)

            else:
                next_frame.append(line)

        return next_frame
                        

                
    def next_step(self):
        self.__changed = []
        marker_p1 = -1
        marker_p2 = -10
        self.losuj(marker_p1, marker_p2)         
        self.__frame = self.new_cell_state(marker_p1, marker_p2)

        self.cykles_counter += 1
            
        return self.__changed

