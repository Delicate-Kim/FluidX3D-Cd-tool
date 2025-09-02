import matplotlib.pyplot as plt

Re = [4054054, 1000000, 405405, 100000, 40541, 10000, 4054, 1000,
      405, 100, 10, 5, 2, 1, 0.5, 0.2, 0.1]
Cd = [0.891019668, 0.8691525, 0.886173999, 0.852312704, 0.863884513,
      0.858254079, 0.864659311, 0.985548121, 1.163738217, 2.041,
      9.25, 15.47377815, 28.786, 43.054, 59.96199856, 77.892, 87.301]

fig, ax = plt.subplots()

ax.set_xscale('log')
ax.set_yscale('log')
ax.scatter(Re, Cd, marker='o') 

ax.set_xlim(1e-1, 1e6) 
ax.set_ylim(6e-2, 4e2)

ax.set_xlabel("Re")
ax.set_ylabel("C_d")
ax.set_title("Drag Coefficient vs Reynolds Number (markers only, log–log)")
ax.grid(True, which='both', ls='--', alpha=0.5)

plt.show()
