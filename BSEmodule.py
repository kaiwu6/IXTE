''' BSE_module'''
from scipy import linalg as LA
import matplotlib as mpl
import matplotlib.pyplot as plt
from const import *
mpl.use('AGG', warn=False)

IND_MAP = {'1': 1, '4': 2, '5': 3}

class SingleLayer():
  '''Class single layer structure'''
  def __init__(self, data_dir, seed_name, _k_grid):
    self.data_dir = data_dir
    self.seed_name = seed_name
    self.k_grid = _k_grid

  def load_data(self):
    '''load data'''
    dim = 5000
    tdim = 0
    band_data = np.zeros((5, 5, dim, 3))
    f = open(self.data_dir + self.seed_name + '_hr.dat', "r")
    f.readline()
    f.readline()
    line = f.readline()
    for i in range(int(int(line.split()[0]) / 15) + 1):
      f.readline()
    while True:
      line = f.readline()
      if line:
        nx, ny, nz, oa, ob, val, = line.split()
        if dist(int(nx), int(ny), 3) and abs(float(val)) > 0.0001 \
              and int(oa) in [1, 4, 5]:
          oa = IND_MAP[oa]
          ob = IND_MAP[ob]
          band_data[oa - 1, ob - 1, tdim, 0] = float(nx)
          band_data[oa - 1, ob - 1, tdim, 1] = float(ny)
          band_data[oa - 1, ob - 1, tdim, 2] = float(val)
          tdim = tdim + 1
      else:
        break
    f.close()
    return tdim, band_data

  def wannier_bands(self):
    tdim, dats = self.load_data()
    # Extract electron and hole band
    eigh = []
    eige = []
    eigd = []
    unitary = []
    unitaryT = []
    for ai in range(NT):
      ka, kb = K_GRID[ai]
      tv = np.zeros((3, 3), dtype=complex)
      for b1 in range(3):
        for b2 in range(3):
          for i in range(tdim):
            ia = dats[b1, b2, i, 0]
            ib = dats[b1, b2, i, 1]
            val = dats[b1, b2, i, 2]
            # ia, ib, val = dats[b1, b2, i]
            tv[b1, b2] += val * np.exp((ia * ka + ib * kb) * 2. * pi * 1j)
      eigs, eigfun = LA.eigh(tv)
      idx = eigs.argsort()[::-1]
      eigs = eigs[idx]
      eigfun = eigfun[:, idx]
      eigh.append(eigs[2].real)
      eige.append(eigs[1].real)
      eigd.append(eigs[0].real)
      unitary.append(eigfun)
      unitaryT.append(np.matrix.getH(eigfun))
      # for Coulomb interaction
    return eigh, eige, eigd, unitary, unitaryT


def dist(_dx, _dy, _n):
  ''' for limited binding distance from first principel calculation'''
  a_dir = np.array([2.0, 0.0])
  b_dir = np.array([-1.0, 1.732])
  r_vec = _dx * a_dir + _dy * b_dir
  tol = 0.1
  return np.sqrt(r_vec[0] ** 2 + r_vec[1] ** 2) <= 2 * _n + tol


def bs_matrix_gen(q_in_1d, xab, eigh, eige, eigd):
  '''generate matrix for each Q of IXs'''
  qa = q_in_1d % KNUMB
  qb = int(q_in_1d / KNUMB)
  dgap = np.min(eige) - np.min(eigh)
  hamilton = np.zeros((2 * NT, 2 * NT), dtype=complex)
  for k1 in range(NT):
    for k2 in range(NT):
      k1q = (k1 + qa) % KNUMB + ((int(k1 / KNUMB) + qb) % KNUMB) * KNUMB
      k2q = (k2 + qa) % KNUMB + ((int(k2 / KNUMB) + qb) % KNUMB) * KNUMB
      dime = k1q + k2q * NT
      dimh = k2 + k1 * NT
      for A2 in range(2):
        for A1 in range(2):
          dim1 = (1 - A1) * NT + k1
          dim2 = (1 - A2) * NT + k2
          hamilton[dim1][dim2] = - V * xab[dime][A1][A2] * xab[dimh][2][2]
  for kb in range(KNUMB):
    for ka in range(KNUMB):
      dim = ka + kb * KNUMB
      dimq = (ka + qa) % KNUMB + ((kb + qb) % KNUMB) * KNUMB
      hamilton[dim][dim] += eige[dimq] - eigh[dim] - dgap
      hamilton[dim + NT][dim + NT] += eigd[dimq] - eigh[dim] - dgap
  return hamilton

def bs_solution(hamiltonian):
  '''solve the function, out put all eigein state and eigen function'''
  eigen_energy, eigen_state, = LA.lapack.zheev(hamiltonian)
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
  x_ab = [np.dot(unitary_t[kk1], unitary[kk2]) for kk1 in range(NT) for kk2 in range(NT)]
  ham = list(map(lambda x: bs_matrix_gen(x, x_ab, eigh, eige, eigd), range(NT)))
  print('diagonalize done')
  ex_results = list(map(bs_solution, ham))
  return ex_results

def fourier_trans(xfun_1):
  '''Fourier transformtion for the wave function'''
  exband = xfun_1[0][0:NT]
  xband = exband.reshape((KNUMB, KNUMB))
  # print('xband', xband)
  x_realspace = np.zeros((KNUMB, KNUMB), dtype=complex)
  for k in range(KNUMB):
    for l in range(KNUMB):
      x_realspace[k][l] = 0
      for i in range(KNUMB):
        for j in range(KNUMB):
          x_realspace[k][l] += xband[i][j] * \
          np.exp((k * i / KNUMB + l * j / KNUMB) * 2 * pi * 1j)
  # print(x_realspace)  # imagine part is zero. good
  return x_realspace

class HighSymmetryPlot():
  '''HighSymmetryPlot'''
  def plot_line_gen(self):
    plot_line = [i for i in range(int(KNUMB / 2))]  # Gamma to M
    plot_line_2 = list(map(lambda x: 2 * KNUMB * x + (int(KNUMB /2) - x), range(int(KNUMB / 6))))  # M to Kbb
    plot_line_3 = list(map(lambda y: (KNUMB + 1) * (int(KNUMB / 3) - y), range(int(KNUMB / 3))))  # K to Gamma
    plot_line.extend(plot_line_2)
    plot_line.extend(plot_line_3)
    return plot_line

  def plot_line_gen_long(self):
    plot_line = [i for i in range(int(KNUMB / 2))]  # Gamma to M
    plot_line_2 = list(map(lambda x: 2 * KNUMB * x + (int(KNUMB/2) - x), range(int(KNUMB / 6))))  # M to K
    plot_line_3 = list(map(lambda y: (KNUMB + 1) * (int(KNUMB / 3) - y), range(int(KNUMB / 3)+1)))  # K to Gamma
    plot_line_4 = list(map(lambda z: (KNUMB + 1) * (KNUMB - z - 1), range(int(KNUMB/3))))  # Gamma to K'
    plot_line_5 = list(map(lambda x: KNUMB * (KNUMB * 2/3 - x - 1) + KNUMB * 2/3 - 4 * (x + 1),
                           range(int(KNUMB / 6))))  # K' to M
    plot_line.extend(plot_line_2)
    plot_line.extend(plot_line_3)
    plot_line.extend(plot_line_4)
    plot_line.extend(plot_line_5)
    return plot_line

  def __init__(self, number_points, data):
    self.number_points = number_points
    self.data = data

  def plot_one_band(self, ax, xeng_1, xeng_2):
    plot_line = self.plot_line_gen()
    num_points = len(plot_line)
    ans = 0
    for ai in range(num_points):
        kxlist = ai
        dim = plot_line[ai]
        ax.scatter(kxlist, xeng_1[dim], s=20.0, facecolors='none',
                   edgecolors='r')
        ax.scatter(kxlist, xeng_2[dim], s=20.0, facecolors='b',
                   edgecolor='none')
    y_max = np.amax([xeng_1, xeng_2]) * 1.05
    ax.set_ylim((0, y_max))
    ax.set_xlim((0, num_points-1))
    plt.show()
    return ans

  def plot_one_band_long(self, ax, xeng_1, xeng_2):
    plot_line = self.plot_line_gen_long()
    num_points = len(plot_line)
    ans = 0
    for ai in range(num_points):
      kxlist = ai
      dim = int(plot_line[ai])
      ax.scatter(kxlist, xeng_1[dim], s=40, facecolors='none',
                   edgecolors='r')
      ax.scatter(kxlist, xeng_2[dim], s=40, facecolors='b',
                   edgecolor='none')
    y_max = np.amax([xeng_1, xeng_2]) * 1.05
    ax.set_ylim((0, y_max))
    ax.set_xlim((0, num_points-1))
    plt.show()
    return ans


def tight_binding(xeng):
  xeng_2d = xeng.reshape((KNUMB, KNUMB))
  hopping = []
  for b in range(3):
    for a in range(b + 1):
      tsum = 0.0
      for j in range(KNUMB):
        for i in range(KNUMB):
          tsum += xeng_2d[i][j] * np.exp((a * i + b * j) * 2. * pi * 1j / KNUMB)
      tsum /= NT
      hopping.append(tsum.real)
  return hopping

def get_band(file_dir):
  return 0

# main program
def part_test():
  '''main'''
  dir_mos2 = '/Users/wk/Dropbox/2015rp/TMDC_Exciton/Tight-binding Model/MoS2/'
  seedname = 'MoS2'
  layer_band = SingleLayer(dir_mos2, seedname, K_GRID)
  eigh, eige, eigd, unitary, unitaryT = layer_band.wannier_bands()
  ex_results = xband(eigh, eige, eigd, unitary, unitaryT)
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
  for i in range(len(hopping)):
    print(hopping[i], file=fo1)
  fo1.close()
  print('step2')

  # tight-binding parameter for exciton2
  hopping = tight_binding(ex_eng_2)
  fo2 = open("t2.dat", "w")
  for i in range(len(hopping)):
    print(hopping[i], file=fo2)
  fo2.close()

  print('start plot')
  fig1 = plt.figure(figsize=(10, 6))
  ax = fig1.add_subplot(111)
  new_plot = HighSymmetryPlot(NT, 1)
  new_plot.plot_one_band_long(ax, ex_eng_1, ex_eng_2)
  plt.savefig('band.png')
  print('end plot')
  '''
  fig2 = plt.figure(figsize=(10, 6))
  ax = fig2.add_subplot(111)
  pcolor(Fouriertrans(xfun))
  '''

if __name__ == "__main__":
    part_test()
