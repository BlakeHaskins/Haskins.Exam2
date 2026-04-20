# region imports
import numpy as np
import math
from scipy.optimize import fsolve
import random as rnd
# endregion

# region class definitions
class UC:
    """
    Unit conversion helper class.
    """

    ft_to_m = 1 / 3.28084
    ft2_to_m2 = ft_to_m ** 2
    ft3_to_m3 = ft_to_m ** 3
    ft3_to_L = ft3_to_m3 * 1000
    L_to_ft3 = 1 / ft3_to_L
    in_to_m = ft_to_m / 12
    m_to_in = 1 / in_to_m
    in2_to_m2 = in_to_m ** 2
    m2_to_in2 = 1 / in2_to_m2
    g_SI = 9.80665
    g_EN = 32.174
    gc_EN = 32.174
    gc_SI = 1.0
    lbf_to_kg = 1 / 2.20462
    lbf_to_N = lbf_to_kg * g_SI
    pa_to_psi = (1 / lbf_to_N) * in2_to_m2

    @classmethod
    def viscosityEnglishToSI(cls, mu, toSI=True):
        """
        Convert viscosity between lb*s/ft^2 and Pa*s.
        """
        # ChatGPT helped write this function.
        cf = (1 / cls.ft2_to_m2) * cls.lbf_to_kg * cls.g_SI
        return mu * cf if toSI else mu / cf

    @classmethod
    def densityEnglishToSI(cls, rho, toSI=True):
        """
        Convert density/specific weight between lb/ft^3 and kg/m^3.
        """
        # ChatGPT helped write this function.
        cf = cls.lbf_to_kg / cls.ft3_to_m3
        return rho * cf if toSI else rho / cf

    @classmethod
    def head_to_pressure(cls, h, rho, SI=True):
        """
        Convert head of fluid to pressure.
        """
        if SI:
            cf = rho * cls.g_SI / cls.gc_SI
            return h * cf
        cf = rho * cls.g_EN / cls.gc_EN * (1 / 12) ** 2
        return h * cf

    @classmethod
    def m_to_psi(cls, h, rho):
        """
        Convert meters of fluid head to psi.
        """
        return cls.head_to_pressure(h, rho) * cls.pa_to_psi

    @classmethod
    def psi_to_m(cls, p, rho):
        """
        Convert psi to meters of fluid head.
        """
        pa = p / cls.pa_to_psi
        return pa / (rho * cls.g_SI)


class Fluid:
    def __init__(self, mu=0.00089, rho=1000, SI=True):
        """
        Store fluid properties internally in SI units.

        Parameters
        ----------
        mu : float
            Dynamic viscosity.
        rho : float
            Density or specific weight input as used by this assignment.
        SI : bool
            True if inputs are already SI.
        """
        # ChatGPT helped write this function.
        self.mu = mu if SI else UC.viscosityEnglishToSI(mu)
        self.rho = rho if SI else UC.densityEnglishToSI(rho)
        self.nu = self.mu / self.rho


class Node:
    def __init__(self, Name='a', Pipes=None, ExtFlow=0):
        """
        A node in the pipe network.

        Parameters
        ----------
        Name : str
            Node label.
        Pipes : list
            Connected Pipe objects.
        ExtFlow : float
            External flow into the node in L/s.
        """
        self.name = Name
        self.pipes = [] if Pipes is None else Pipes
        self.extFlow = ExtFlow
        self.QNet = 0.0
        self.P = 0.0
        self.oCalculated = False

    def getNetFlowRate(self):
        """
        Compute net flow into the node in L/s.
        """
        qtot = self.extFlow
        for p in self.pipes:
            qtot += p.getFlowIntoNode(self.name)
        self.QNet = qtot
        return self.QNet

    def setExtFlow(self, E, SI=True):
        """
        Set node external flow.

        Parameters
        ----------
        E : float
            External flow value.
        SI : bool
            False means E is in cfs and should be converted to L/s.
        """
        self.extFlow = E if SI else E * UC.ft3_to_L


class Loop:
    def __init__(self, Name='A', Pipes=None):
        """
        Pipe loop listed in traversal order.
        """
        self.name = Name
        self.pipes = [] if Pipes is None else Pipes

    def getLoopHeadLoss(self):
        """
        Compute signed head change around the loop in meters of water.
        """
        deltaP = 0.0
        startNode = self.pipes[0].startNode
        for p in self.pipes:
            deltaP += p.getFlowHeadLoss(startNode)
            startNode = p.endNode if startNode != p.endNode else p.startNode
        return deltaP


class Pipe:
    def __init__(self, Start='A', End='B', L=100, D=200, r=0.00025, fluid=Fluid(), SI=True):
        """
        Pipe object.

        Parameters
        ----------
        Start : str
            One endpoint.
        End : str
            Other endpoint.
        L : float
            Length in m if SI=True, otherwise ft.
        D : float
            Diameter in mm if SI=True, otherwise inches.
        r : float
            Roughness in m if SI=True, otherwise ft.
        fluid : Fluid
            Fluid object.
        SI : bool
            False means convert English inputs to SI internally.
        """
        # ChatGPT helped write this function.
        self.startNode = min(Start.lower(), End.lower())
        self.endNode = max(Start.lower(), End.lower())
        self.length = L if SI else UC.ft_to_m * L
        self.rough = r if SI else UC.ft_to_m * r
        self.fluid = fluid

        self.d = D / 1000.0 if SI else UC.in_to_m * D
        self.relrough = self.rough / self.d
        self.A = math.pi * self.d ** 2 / 4.0

        self.Q = 10.0   # L/s initial guess
        self.vel = self.V()
        self.reynolds = self.Re()
        self.hl = 0.0

    def V(self):
        """
        Average velocity in m/s.
        """
        self.vel = (self.Q / 1000.0) / self.A
        return self.vel

    def Re(self):
        """
        Reynolds number.
        """
        self.reynolds = self.fluid.rho * abs(self.V()) * self.d / self.fluid.mu
        return self.reynolds

    def FrictionFactor(self):
        """
        Darcy friction factor.
        """
        Re = self.Re()
        rr = self.relrough

        def colebrook():
            fn = lambda f: 1 / np.sqrt(f) + 2.0 * np.log10(rr / 3.7 + 2.51 / (Re * np.sqrt(f)))
            result = fsolve(fn, 0.01)
            return result[0]

        def laminar():
            return 64.0 / Re

        if Re >= 4000:
            return colebrook()
        if Re <= 2000:
            return laminar()

        ff_cb = colebrook()
        ff_lam = laminar()
        mean = ff_lam + (Re - 2000.0) * (ff_cb - ff_lam) / 2000.0

        sig_1 = (1 - (Re - 3000) / 1000) * 0.2 * mean
        sig_2 = (1 - (3000 - Re) / 1000) * 0.2 * mean
        sig = sig_1 if Re >= 3000 else sig_2
        return rnd.normalvariate(mean, sig)

    def frictionHeadLoss(self):
        """
        Darcy-Weisbach head loss in meters of water.
        Always positive magnitude.
        """
        g = 9.81
        ff = self.FrictionFactor()
        self.hl = ff * (self.length / self.d) * (abs(self.V()) ** 2) / (2 * g)
        return self.hl

    def getFlowHeadLoss(self, s):
        """
        Signed head loss for loop traversal.
        """
        nTraverse = 1 if s == self.startNode else -1
        nFlow = 1 if self.Q >= 0 else -1
        return nTraverse * nFlow * self.frictionHeadLoss()

    def Name(self):
        """
        Pipe name as start-end.
        """
        return self.startNode + '-' + self.endNode

    def oContainsNode(self, node):
        """
        True if the pipe connects to the given node.
        """
        return self.startNode == node or self.endNode == node

    def printPipeFlowRate(self, SI=True):
        """
        Print pipe flow and Reynolds number.
        """
        q_units = 'L/s' if SI else 'cfs'
        q = self.Q if SI else self.Q * UC.L_to_ft3
        print('The flow in segment {} is {:0.2f} ({}) and Re={:.1f}'.format(
            self.Name(), q, q_units, abs(self.reynolds)))

    def printPipeHeadLoss(self, SI=True):
        """
        Print pipe head loss.
        """
        cfd = 1000 if SI else UC.m_to_in
        unitsd = 'mm' if SI else 'in'
        cfL = 1 if SI else 1 / UC.ft_to_m
        unitsL = 'm' if SI else 'ft'
        cfh = cfd
        units_h = unitsd
        print("head loss in pipe {} (L={:.2f} {}, d={:.2f} {}) is {:.2f} {} of water".format(
            self.Name(), self.length * cfL, unitsL, self.d * cfd, unitsd, self.hl * cfh, units_h))

    def getFlowIntoNode(self, n):
        """
        Signed flow into node n.
        """
        if n == self.startNode:
            return -self.Q
        return self.Q


class PipeNetwork:
    def __init__(self, Pipes=None, Loops=None, Nodes=None, fluid=Fluid()):
        """
        Pipe network object.
        """
        self.loops = [] if Loops is None else Loops
        self.nodes = [] if Nodes is None else Nodes
        self.Fluid = fluid
        self.pipes = [] if Pipes is None else Pipes

    def findFlowRates(self):
        """
        Solve for all pipe flow rates.

        Returns
        -------
        ndarray
            Pipe flow rates in L/s.
        """
        # one unknown per pipe
        N = len(self.pipes)

        # better initial guess based on the example output on the exam page
        # example values are in cfs, so convert to L/s
        Q0 = np.array([3.57, -3.57, 2.57, 1.01, 0.53, 2.04, -1.47,
                       1.46, -3.45, -1.50, -2.97, 6.43, 2.97]) * UC.ft3_to_L

        def fn(q):
            """
            Residual vector for fsolve.
            """
            for n in self.nodes:
                n.P = 0.0
                n.oCalculated = False

            for i in range(len(self.pipes)):
                self.pipes[i].Q = q[i]

            node_res = self.getNodeFlowRates()

            # drop one redundant node equation so total equations = number of pipes
            node_res = node_res[:-1]

            loop_res = self.getLoopHeadLosses()

            return np.array(node_res + loop_res)

        FR = fsolve(fn, Q0)

        for i in range(len(self.pipes)):
            self.pipes[i].Q = FR[i]

        self.getNodeFlowRates()
        self.getLoopHeadLosses()

        return FR

    def getNodeFlowRates(self):
        """
        Net flow into each node.
        """
        return [n.getNetFlowRate() for n in self.nodes]

    def getLoopHeadLosses(self):
        """
        Net head loss for each loop.
        """
        return [l.getLoopHeadLoss() for l in self.loops]

    def getNodePressures(self, knownNodeP, knownNode):
        """
        Compute node pressures from one known reference node.
        """
        for n in self.nodes:
            n.P = 0.0
            n.oCalculated = False

        for l in self.loops:
            startNode = l.pipes[0].startNode
            n = self.getNode(startNode)
            currentP = n.P
            n.oCalculated = True

            for p in l.pipes:
                phl = p.getFlowHeadLoss(startNode)
                currentP -= phl
                startNode = p.endNode if startNode != p.endNode else p.startNode
                n = self.getNode(startNode)
                n.P = currentP

        kn = self.getNode(knownNode)
        deltaP = knownNodeP - kn.P
        for n in self.nodes:
            n.P += deltaP

    def getPipe(self, name):
        """
        Return pipe by name.
        """
        for p in self.pipes:
            if name == p.Name():
                return p

    def getNodePipes(self, node):
        """
        Return all pipes connected to a node.
        """
        out = []
        for p in self.pipes:
            if p.oContainsNode(node):
                out.append(p)
        return out

    def nodeBuilt(self, node):
        """
        Check whether a node already exists.
        """
        for n in self.nodes:
            if n.name == node:
                return True
        return False

    def getNode(self, name):
        """
        Return node by name.
        """
        for n in self.nodes:
            if n.name == name:
                return n

    def buildNodes(self):
        """
        Build nodes automatically from pipes.
        """
        for p in self.pipes:
            if not self.nodeBuilt(p.startNode):
                self.nodes.append(Node(p.startNode, self.getNodePipes(p.startNode)))
            if not self.nodeBuilt(p.endNode):
                self.nodes.append(Node(p.endNode, self.getNodePipes(p.endNode)))

    def printPipeFlowRates(self, SI=True):
        for p in self.pipes:
            p.printPipeFlowRate(SI=SI)

    def printNetNodeFlows(self, SI=True):
        for n in self.nodes:
            Q = n.QNet if SI else n.QNet * UC.L_to_ft3
            units = 'L/S' if SI else 'cfs'
            print('net flow into node {} is {:0.2f} ({})'.format(n.name, Q, units))

    def printLoopHeadLoss(self, SI=True):
        cf = UC.m_to_psi(1, self.pipes[0].fluid.rho)
        units = 'm of water' if SI else 'psi'
        for l in self.loops:
            hl = l.getLoopHeadLoss()
            hl = hl if SI else hl * cf
            print('head loss for loop {} is {:0.2f} ({})'.format(l.name, hl, units))

    def printPipeHeadLoss(self, SI=True):
        for p in self.pipes:
            p.printPipeHeadLoss(SI=SI)

    def printNodePressures(self, SI=True):
        pUnits = 'm of water' if SI else 'psi'
        cf = 1.0 if SI else UC.m_to_psi(1, self.Fluid.rho)
        for n in self.nodes:
            p = n.P * cf
            print('Pressure at node {} = {:0.2f} {}'.format(n.name, p, pUnits))
# endregion


# region function definitions
def main():
    """
    Solve MAE 3403 Exam 2 Problem 1 pipe network.

    Conventions:
    1. Positive pipe flow is from lower-letter node to higher-letter node.
    2. Mass is conserved at each node.
    3. Net head loss around each loop is zero.
    4. Node h pressure is known to be 80 psi.
    """
    SIUnits = False

    # room-temperature water from the problem statement
    water = Fluid(mu=20.50e-6, rho=62.3, SI=SIUnits)

    # roughness values from the problem statement
    r_CI = 0.00085   # cast iron for 12 in and 16 in pipes
    r_CN = 0.003     # concrete for 18 in and 24 in pipes

    PN = PipeNetwork()
    PN.Fluid = water

    # pipe definitions from the figure
    PN.pipes.append(Pipe('a', 'b', 1000, 18, r_CN, water, SI=SIUnits))
    PN.pipes.append(Pipe('a', 'h', 1600, 24, r_CN, water, SI=SIUnits))
    PN.pipes.append(Pipe('b', 'c', 500, 18, r_CN, water, SI=SIUnits))
    PN.pipes.append(Pipe('b', 'e', 800, 16, r_CI, water, SI=SIUnits))
    PN.pipes.append(Pipe('c', 'd', 500, 18, r_CN, water, SI=SIUnits))
    PN.pipes.append(Pipe('c', 'f', 800, 16, r_CI, water, SI=SIUnits))
    PN.pipes.append(Pipe('d', 'g', 800, 16, r_CI, water, SI=SIUnits))
    PN.pipes.append(Pipe('e', 'f', 500, 12, r_CI, water, SI=SIUnits))
    PN.pipes.append(Pipe('e', 'i', 800, 18, r_CN, water, SI=SIUnits))
    PN.pipes.append(Pipe('f', 'g', 500, 12, r_CI, water, SI=SIUnits))
    PN.pipes.append(Pipe('g', 'j', 800, 18, r_CN, water, SI=SIUnits))
    PN.pipes.append(Pipe('h', 'i', 1000, 24, r_CN, water, SI=SIUnits))
    PN.pipes.append(Pipe('i', 'j', 1000, 24, r_CN, water, SI=SIUnits))

    PN.buildNodes()

    # external flows from the exam figure
    PN.getNode('h').setExtFlow(10, SI=SIUnits)
    PN.getNode('e').setExtFlow(-3, SI=SIUnits)
    PN.getNode('f').setExtFlow(-5, SI=SIUnits)
    PN.getNode('d').setExtFlow(-2, SI=SIUnits)

    # loops in continuous traversal order
    PN.loops.append(Loop('A', [
        PN.getPipe('a-b'),
        PN.getPipe('b-e'),
        PN.getPipe('e-i'),
        PN.getPipe('h-i'),
        PN.getPipe('a-h')
    ]))

    PN.loops.append(Loop('B', [
        PN.getPipe('b-c'),
        PN.getPipe('c-f'),
        PN.getPipe('e-f'),
        PN.getPipe('b-e')
    ]))

    PN.loops.append(Loop('C', [
        PN.getPipe('c-d'),
        PN.getPipe('d-g'),
        PN.getPipe('f-g'),
        PN.getPipe('c-f')
    ]))

    PN.loops.append(Loop('D', [
        PN.getPipe('e-i'),
        PN.getPipe('i-j'),
        PN.getPipe('g-j'),
        PN.getPipe('f-g'),
        PN.getPipe('e-f')
    ]))

    PN.findFlowRates()

    knownP = UC.psi_to_m(80, water.rho)
    PN.getNodePressures(knownNode='h', knownNodeP=knownP)

    PN.printPipeFlowRates(SI=SIUnits)
    print()
    print('Check node flows:')
    PN.printNetNodeFlows(SI=SIUnits)
    print()
    print('Check loop head loss:')
    PN.printLoopHeadLoss(SI=SIUnits)
    print()
    PN.printPipeHeadLoss(SI=SIUnits)
    print()
    PN.printNodePressures(SI=SIUnits)
# endregion


# region function calls
if __name__ == "__main__":
    main()
# endregion