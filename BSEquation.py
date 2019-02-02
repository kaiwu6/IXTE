''' BSEquation
calculate IX bands from BS equations:
  input: directory and seedname of the single-layer files
  output:  first and second energy bands of IX (two array)
1. Load the result from single layer and energy spectrums
2. Build the BS matrix from interaction and energy spectrums
3. Diagonalize the BS matrix to get first and second bands
4. Band plot along high symmetry points.
'''
import numpy as np
from scipy import linalg as LA
#import matplotlib as mpl
import matplotlib.pyplot as plt
from const import pi, KNUMB, NT, K_ARRAY, V
#mpl.use('AGG', warn=False)

IND_MAP = {'1': 1, '4': 2, '5': 3}
class SingleLayer():
  '''Class single layer structure'''
  def __init__(self, data_dir, seed_name):
    self.data_dir = data_dir
    self.seed_name = seed_name

  def load_data(self):
    '''load data'''
    dim = 5000
    tdim = 0
    band_data = np.zeros((5, 5, dim, 3))
    f_hr = open(self.data_dir + self.seed_name + '_hr.dat', "r")
    f_hr.readline()
    f_hr.readline()
    line = f_hr.readline()
    for _ in range(int(int(line.split()[0]) / 15) + 1):
      f_hr.readline()
    while True:
      line = f_hr.readline()
      if line:
        n_x, n_y, _, o_a, o_b, val, _ = line.split()
        if dist(int(n_x), int(n_y), 3) and abs(float(val)) > 0.0001 \
              and int(o_a) in [1, 4, 5]:
          o_a = IND_MAP[o_a]
          o_b = IND_MAP[o_b]
          band_data[o_a - 1, o_b - 1, tdim, 0] = float(n_x)
          band_data[o_a - 1, o_b - 1, tdim, 1] = float(n_y)
          band_data[o_a - 1, o_b - 1, tdim, 2] = float(val)
          tdim = tdim + 1
      else:
        break
    f_hr.close()
    return tdim, band_data

  def wannier_bands(self):
    '''Diagonalize and extract electron(2) and hole(1) bands'''
    tdim, dats = self.load_data()
    eigh = []
    eige = []
    eigd = []
    unitary = []
    unitary_trans = []
    for a_i in range(NT):
      k_a, k_b = K_ARRAY[a_i]
      t_v = np.zeros((3, 3), dtype=complex)
      for b_1 in range(3):
        for b_2 in range(3):
          for i in range(tdim):
            i_a = dats[b_1, b_2, i, 0]
            i_b = dats[b_1, b_2, i, 1]
            val = dats[b_1, b_2, i, 2]
            t_v[b_1, b_2] += val * np.exp((
                i_a * k_a + i_b * k_b) * 1j)
      eigs, eigfun = LA.eigh(t_v)
      idx = eigs.argsort()[::-1]
      eigs = eigs[idx]
      eigfun = eigfun[:, idx]
      eigh.append(eigs[2].real)
      eige.append(eigs[1].real)
      eigd.append(eigs[0].real)
      #Eigen function is for the Coulomb interaction
      unitary.append(eigfun)
      unitary_trans.append(np.matrix.getH(eigfun))
    return eigh, eige, eigd, unitary, unitary_trans

def dist(_dx, _dy, _n1):
  ''' for limied binding distance from first principel calculation'''
  a_dir = np.array([2.0, 0.0])
  b_dir = np.array([-1.0, 1.732])
  r_vec = _dx * a_dir + _dy * b_dir
  tol = 0.1
  return np.sqrt(r_vec[0] ** 2 + r_vec[1] ** 2) <= 2 * _n1 + tol

def bs_matrix_gen(q_in_1d, xab, eigh, eige, eigd):
  '''generate matrix for each Q of IXs'''
  q_a = q_in_1d % KNUMB
  q_b = int(q_in_1d / KNUMB)
  gap = np.min(eige) - np.min(eigh)
  hamilton = np.zeros((2 * NT, 2 * NT), dtype=complex)
  for k_1 in range(NT):
    for k_2 in range(NT):
      k1q = (k_1 + q_a) % KNUMB + ((int(k_1 / KNUMB) + q_b) % KNUMB) * KNUMB
      k2q = (k_2 + q_a) % KNUMB + ((int(k_2 / KNUMB) + q_b) % KNUMB) * KNUMB
      dime = k1q + k2q * NT
      dimh = k_2 + k_1 * NT
      for electron in range(2):
        for hole in range(2):
          dim1 = (1 - hole) * NT + k_1
          dim2 = (1 - electron) * NT + k_2
          hamilton[dim1][dim2] -= V*xab[dime][hole][electron]*xab[dimh][2][2]
  for k_b in range(KNUMB):
    for k_a in range(KNUMB):
      dim = k_a + k_b * KNUMB
      dimq = (k_a + q_a) % KNUMB + ((k_b + q_b) % KNUMB) * KNUMB
      hamilton[dim][dim] += eige[dimq] - eigh[dim] - gap
      hamilton[dim + NT][dim + NT] += eigd[dimq] - eigh[dim] - gap
  return hamilton

def bs_solution(hamiltonian):
  '''solve the function, out put all eigein state and eigen function'''
  eigen_energy, eigen_state, _ = LA.lapack.zheev(hamiltonian)
  indx = eigen_energy.argsort()
  eigen_energy = eigen_energy[indx]
  eigen_state = eigen_state[:, indx]
  x_eng_1 = eigen_energy[0].real
  x_fun_1 = eigen_state[:, 0]
  x_eng_2 = eigen_energy[1].real
  x_fun_2 = eigen_state[:, 1]
  return x_eng_1, x_fun_1, x_eng_2, x_fun_2

def xband(eigh, eige, eigd, unitary, unitary_t):
  '''indirect exction strcuture'''
  x_ab = [np.dot(unitary_t[kk1], unitary[kk2])
          for kk1 in range(NT) for kk2 in range(NT)]
  ham = list(map(lambda x: bs_matrix_gen(x, x_ab, eigh, eige, eigd), range(NT)))
  print('diagonalize done')
  ex_results = list(map(bs_solution, ham))
  return ex_results


def fourier_trans(xfun_1):
  '''Fourier transformtion for the wave function'''
  exband = xfun_1[0][0:NT]
  ix_band = exband.reshape((KNUMB, KNUMB))
  # print('xband', xband)
  x_realspace = np.zeros((KNUMB, KNUMB), dtype=complex)
  for k_x in range(KNUMB):
    for k_y in range(KNUMB):
      x_realspace[k_x][k_y] = 0
      for i in range(KNUMB):
        for j in range(KNUMB):
          x_realspace[k_x][k_y] += ix_band[i][j] * \
          np.exp((k_x * i / KNUMB + k_y * j / KNUMB) * 2 * pi * 1j)
  # print(x_realspace)  # imagine part is zero. good
  return x_realspace

class HighSymmetryPlot():
  '''HighSymmetryPlot'''
  def plot_line_gen(self):
    '''plot line'''
    plot_line = [i for i in range(int(KNUMB / 2))]  # Gamma to M
    plot_line_2 = list(map(lambda x: 2*KNUMB*x + (int(KNUMB/2)-x),
                           range(int(KNUMB/6))))  # M to Kbb
    plot_line_3 = list(map(lambda y: (KNUMB+1) * (int(KNUMB/3)-y),
                           range(int(KNUMB/3))))  # K to Gamma
    plot_line.extend(plot_line_2)
    plot_line.extend(plot_line_3)
    return plot_line

  def plot_line_gen_long(self):
    ''' long plot '''
    plot_line = [i for i in range(int(KNUMB / 2))]  # Gamma to M
    plot_line_2 = list(map(lambda x: 2 * KNUMB * x + (int(KNUMB/2) - x),
                           range(int(KNUMB / 6))))  # M to K
    plot_line_3 = list(map(lambda y: (KNUMB + 1) * (int(KNUMB / 3) - y),
                           range(int(KNUMB / 3)+1)))  # K to Gamma
    plot_line_4 = list(map(lambda z: (KNUMB + 1) * (KNUMB - z - 1),
                           range(int(KNUMB/3))))  #Gamma to K'
    plot_line_5 = list(map(lambda x: KNUMB * (KNUMB * 2/3 - x - 1)
                           +KNUMB * 2/3 - 4 * (x + 1), range(int(KNUMB / 6))))
                           #K' to M
    plot_line.extend(plot_line_2)
    plot_line.extend(plot_line_3)
    plot_line.extend(plot_line_4)
    plot_line.extend(plot_line_5)
    return plot_line

  def __init__(self, number_points, data):
    self.number_points = number_points
    self.data = data

  def plot_one_band(self, a_x, xeng_1, xeng_2):
    ''' band plot one '''
    plot_line = self.plot_line_gen()
    num_points = len(plot_line)
    ans = 0
    for a_i in range(num_points):
      kxlist = a_i
      dim = plot_line[a_i]
      a_x.scatter(kxlist, xeng_1[dim], s=20.0, facecolors='none',
                  edgecolors='r')
      a_x.scatter(kxlist, xeng_2[dim], s=20.0, facecolors='b',
                  edgecolor='none')
    y_max = np.amax([xeng_1, xeng_2]) * 1.05
    a_x.set_ylim((0, y_max))
    a_x.set_xlim((0, num_points-1))
    plt.show()
    return ans

  def plot_one_band_long(self, a_x, xeng_1, xeng_2):
    '''band plot long'''
    plot_line = self.plot_line_gen_long()
    num_points = len(plot_line)
    ans = 0
    for a_i in range(num_points):
      kxlist = a_i
      dim = int(plot_line[a_i])
      a_x.scatter(kxlist, xeng_1[dim], s=40, facecolors='none',
                  edgecolors='r')
      a_x.scatter(kxlist, xeng_2[dim], s=40, facecolors='b',
                  edgecolor='none')
    y_max = np.amax([xeng_1, xeng_2]) * 1.05
    a_x.set_ylim((0, y_max))
    a_x.set_xlim((0, num_points-1))
    plt.show()
    return ans

def tight_binding(xeng):
  ''' tight_banding test '''
  xeng_2d = xeng.reshape((KNUMB, KNUMB))
  hopping = []
  for t_b in range(3):
    for t_a in range(t_b + 1):
      tsum = 0.0
      for j in range(KNUMB):
        for i in range(KNUMB):
          tsum += xeng_2d[i][j] * np.exp((t_a*i+t_b*j)*2.*pi*1j/KNUMB)
      tsum /= NT
      hopping.append(tsum.real)
  return hopping

def compute_ixbands(file_dir, seed_name):
  '''callable function for the final output'''
  layer_band = SingleLayer(file_dir, seed_name)
  eigh, eige, eigd, unitary, unitary_trans = layer_band.wannier_bands()
  ix_results = xband(eigh, eige, eigd, unitary, unitary_trans)
  ix_band_1 = list(map(lambda x: ix_results[x][0], range(NT)))
  ix_band_2 = list(map(lambda x: ix_results[x][2], range(NT)))
  ground_state = np.amin([ix_band_1, ix_band_2])
  #print(ground_state)
  #print(np.amin(ix_band_1))
  #print(np.amin(ix_band_2))
  ix_band_1 -= ground_state
  ix_band_2 -= ground_state
  print('band calulation completed: ground state energy = 0')
  return ix_band_1, ix_band_2

# main program
def unit_test():
  '''main'''
  dir_mos2 = '/Users/wk/Dropbox/2015rp/TMDC_Exciton/Tight-binding Model/MoS2/'
  seedname = 'MoS2'
  layer_band = SingleLayer(dir_mos2, seedname)
  eigh, eige, eigd, unitary, unitary_trans = layer_band.wannier_bands()
  ex_results = xband(eigh, eige, eigd, unitary, unitary_trans)
  ex_eng_1 = list(map(lambda x: ex_results[x][0], range(NT)))
  #ex_fun_1 = list(map(lambda x: ex_results[x][1], range(NT)))
  ex_eng_2 = list(map(lambda x: ex_results[x][2], range(NT)))
  #ex_fun_2 = list(map(lambda x: ex_results[x][3], range(NT)))
  xengmin = np.amin(ex_eng_1)
  # print('ex_min', xengmin)
  ex_eng_1 -= xengmin
  ex_eng_2 -= xengmin
  # print(ex_eng_1, ex_eng_2)
  print("step 1: get all band")

  # tight-binding parameter for excitons
  hopping = tight_binding(ex_eng_1)
  fo1 = open("t1.dat", "w")
  for hop in hopping:
    print(hop, file=fo1)
  fo1.close()
  print('step2')

  # tight-binding parameter for exciton2
  hopping = tight_binding(ex_eng_2)
  fo2 = open("t2.dat", "w")
  for hop in hopping:
    print(hop, file=fo2)
  fo2.close()

  print('start plot')
  fig1 = plt.figure(figsize=(10, 6))
  a_x = fig1.add_subplot(111)
  new_plot = HighSymmetryPlot(NT, 1)
  new_plot.plot_one_band_long(a_x, ex_eng_1, ex_eng_2)
  plt.savefig('band.png')
  print('end plot')

if __name__ == "__main__":
  unit_test()
