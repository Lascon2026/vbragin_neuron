from neuron import h, gui
import matplotlib.pyplot as plt

class HHCell: 
    """Two-section cell: A soma with active channels and
    a dendrite with passive properties."""
    def __init__(self):
        self.synlist = []
        self.nclist = []
        self.nslist = []

        self.create_sections()
        self.build_topology()
        self.define_geometry()
        self.define_biophysics()

    def create_sections(self):
        """Create the sections of the cell."""
        self.soma = h.Section(name='soma')
        self.dend = h.Section(name='dend')
    
    def build_topology(self):
        """Connect the sections of the cell"""
        self.dend.connect(self.soma(1))
    
    def define_geometry(self):
        """Set the 3D geometry of the cell."""
        self.soma.L = self.soma.diam = 12.6157 # microns
        self.dend.L = 200                      # microns
        self.dend.diam = 1                     # microns
        self.dend.nseg = 9
    
    def define_biophysics(self):
        """Assign the membrane properties across the cell."""
        for sec in [self.soma, self.dend]: # 
            sec.Ra = 100    # Axial resistance in Ohm * cm
            sec.cm = 1      # Membrane capacitance in micro Farads / cm^2
        
        # Insert active Hodgkin-Huxley current in the soma
        self.soma.insert('hh')
        self.soma.gnabar_hh = 0.12  # Sodium conductance in S/cm2
        self.soma.gkbar_hh = 0.036  # Potassium conductance in S/cm2
        self.soma.gl_hh = 0.0003    # Leak conductance in S/cm2
        self.soma.el_hh = -54.3     # Reversal potential in mV
        # Insert passive current in the dendrite
        self.dend.insert('pas')
        self.dend.g_pas = 0.001  # Passive conductance in S/cm2
        self.dend.e_pas = -65    # Leak reversal potential mV

    def add_current_stim(self, delay):
        self.stim = h.IClamp(self.dend(1.0))
        self.stim.amp = 0.3  # input current in nA
        self.stim.delay = delay  # turn on after this time in ms
        self.stim.dur = 1  # duration of 1 ms

    def set_recording(self):
        """Set soma, dendrite, and time recording vectors on the cell. """
        self.soma_v_vec = h.Vector()   # Membrane potential vector at soma
        self.dend_v_vec = h.Vector()   # Membrane potential vector at dendrite
        self.t_vec = h.Vector()        # Time stamp vector
        self.soma_v_vec.record(self.soma(0.5)._ref_v)
        self.dend_v_vec.record(self.dend(0.5)._ref_v)
        self.t_vec.record(h._ref_t)

    def plot_voltage(self, title='Cell voltage', ylim=None, show=True):
        """Plot the recorded traces"""
        fig = plt.figure(figsize=(8,4)) # Default figsize is (8,6)
        plt.plot(self.t_vec, self.soma_v_vec, color='black', label='soma(0.5)')
        plt.plot(self.t_vec, self.dend_v_vec, color='red', label='dend(0.5)')
        plt.legend()
        plt.xlabel('time (ms)')
        plt.ylabel('mV')
        plt.ylim(ylim)
        plt.title(title)
        if show:
            plt.show()
        return fig

    def create_synapse(self, loc=0.5, tau=2, e=0):
        syn = h.ExpSyn(self.dend(loc))
        syn.tau = tau
        syn.e = e
        self.synlist.append(syn)
    
    def connect2pre(self, preCell, synid=0, delay=2, weight=1):
        nc = h.NetCon(preCell.soma(0.5)._ref_v, self.synlist[synid], sec = preCell.soma)
        nc.delay = delay
        nc.weight[0] = weight
        self.nclist.append(nc)

    def set_position(self, x, y):
        self.x = x
        self.y = y


import numpy as np
np.random.seed(111) # set fixed "random seed" to ensure reproducibility across the trials

class Pop:

    size = (100, 100) # width, height

    def __init__(self, numCells, xNormRange=[0.0, 1.0], yNormRange=[0.0, 1.0]):
        self.spkt = h.Vector()   # Spike time of all cells
        self.spkid = h.Vector()  # cell ids of spike times
        self.create_cells(numCells, xNormRange, yNormRange)  # call method to create cells
 
    def create_cells(self, numCells, xNormRange, yNormRange):
        """ Create cells in the network """
        self.cells = []

        xMin = self.size[0] * xNormRange[0]
        xMax = self.size[0] * xNormRange[1]
        xPoss = np.random.uniform(xMin, xMax, numCells)

        yMin = self.size[1] * yNormRange[0]
        yMax = self.size[1] * yNormRange[1]
        yPoss = np.random.uniform(yMin, yMax, numCells)

        for i in range(numCells):  # for each cell
            cell = HHCell()  # create cell object
            cell.set_recording()  # set up voltage recording
            cell.create_synapse()

            cell.set_position(xPoss[i], yPoss[i])  # set x,y position
            nc = h.NetCon(cell.soma(0.5)._ref_v, None, sec=cell.soma)  # create netcon to record spikes from cell
            nc.record(self.spkt, self.spkid, i) 
            self.cells.append(cell)  # add cell to list of cells in network
            print('Created cell ' + str(i))

    def plot_net(self, show=True):
        """ Plot position of cells and conns """

        plt.xlim(0, self.size[0])
        plt.ylim(0, self.size[1])
        posX = [cell.x for cell in self.cells]  # get all x positions
        posY = [cell.y for cell in self.cells]  # get all y positions
        plt.scatter(posX, posY, s=40) # plot cell soma positions
        if show:
            plt.show()

    def plot_raster(self, color='blue', show=True):
        """ Plot raster with spikes of all cells """
        plt.figure()
        plt.scatter(list(self.spkt), list(self.spkid), marker= "|", s=100, c=color)
        plt.xlabel('time (ms)')
        plt.ylabel('cell id')
        plt.title('Network raster')
        if show:
            plt.show()

popExc = Pop(100, yNormRange=[0.2, 1])
popInh = Pop(40, yNormRange=[0, 0.33])

# plot cells' locations
popExc.plot_net(show=False)
popInh.plot_net(show=False)

h.tstop = 1e3 # set simulation duration
h.init()
h.run()  # run simulation

# plot raster for both pops
popExc.plot_raster(show=False)
# popInh.plot_raster(show=False) # it's empty so far

# plot voltage of some selected cell
popExc.cells[0].plot_voltage(show=False)

plt.show()
