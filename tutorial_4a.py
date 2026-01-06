from neuron import h, gui
import matplotlib.pyplot as plt

class HHCell: 
    """Two-section cell: A soma with active channels and
    a dendrite with passive properties."""
    def __init__(self):
        self.synlist = []

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
        self.dend.nseg = 10
    
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

    def create_synapse(self, e, loc=0.5, tau=2):
        syn = h.ExpSyn(self.dend(loc))
        syn.tau = tau
        syn.e = e
        self.synlist.append(syn)
    
    def connect2pre(self, preCell, synid=0, delay=2, weight=1):
        nc = h.NetCon(preCell.soma(0.5)._ref_v, self.synlist[synid], sec = preCell.soma)
        nc.delay = delay
        nc.weight[0] = weight

# create cells
preCell1 = HHCell()   # create presyn cell 1
preCell1.add_current_stim(20)  # add stimulation
preCell1.set_recording()  # setup recording

preCell2 = HHCell()  # create presyn cell 2
preCell2.add_current_stim(20)  # add stimulation
preCell2.set_recording()  # setup recording

postCell = HHCell()  # create postsyn cell
postCell.set_recording()  # setup recording

# create synapses

# .....


# connect cells

# .....


h.tstop = 80 # set simulation duration
h.init()
h.run()  # run simulation

preCell1.plot_voltage(show=False, title='preCell1 voltage')  # plot voltage
preCell2.plot_voltage(show=False, title='preCell2 voltage')  # plot voltage
postCell.plot_voltage(show=False, title='postCell voltage', ylim=(-80, 30))  # plot voltage

plt.show()
