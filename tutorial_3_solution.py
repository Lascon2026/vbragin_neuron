from neuron import h, gui
import matplotlib.pyplot as plt

class HHCell: 
    """Two-section cell: A soma with active channels and
    a dendrite with passive properties."""
    def __init__(self):
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


postCell = HHCell()
postCell.set_recording()

preCell1 = HHCell()
preCell1.set_recording()
preCell1.add_current_stim(delay=10)

preCell2 = HHCell()
preCell2.set_recording()
preCell2.add_current_stim(delay=10)


# connect preCell1 -> postCell to generate small EPSP (or spike)
postSyn1 = h.ExpSyn(postCell.dend(0.5))
postSyn1.e = 0 # excitatory synapse
postSyn1.tau = 2

pre1Conn = h.NetCon(preCell1.soma(0.5)._ref_v, postSyn1, sec=preCell1.soma)
# pre1Conn.weight[0] = 0.002 # small EPSP
pre1Conn.weight[0] = 0.005 # triggers spike postsynaptically
pre1Conn.delay = 1


# # connect preCell2 -> postCell to generate another small EPSP
# postSyn2 = h.ExpSyn(postCell.dend(0.5))
# postSyn2.e = 0 # exciatory synapse
# postSyn2.tau = 2

# pre2Conn = h.NetCon(preCell2.soma(0.5)._ref_v, postSyn2, sec=preCell2.soma)
# pre2Conn.weight[0] = 0.002
# pre2Conn.delay = 1


# Alternatively, connect preCell2 -> postCell to generate IPSP
postSyn2 = h.ExpSyn(postCell.dend(0.5))
postSyn2.e = -80 # inhibitory synapse
postSyn2.tau = 2

pre2Conn = h.NetCon(preCell2.soma(0.5)._ref_v, postSyn2, sec=preCell2.soma)
pre2Conn.weight[0] = 0.01
pre2Conn.delay = 1



h.tstop = 30 # set simulation duration
h.init()
h.run()  # run simulation

# preCell1.plot_voltage(show=False, title='preCell1 voltage')  # plot voltage
# preCell2.plot_voltage(show=False, title='preCell2 voltage')  # plot voltage
postCell.plot_voltage(show=False, title='postCell voltage', ylim=(-80, 30))  # plot voltage

plt.show()