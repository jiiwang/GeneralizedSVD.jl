# import matplotlib.pyplot as plt
# import matplotlib.colors as colors
# import numpy as np
#
# # This dictionary defines the colormap
# cdict = {'red':  ((0.0, 0.0, 0.0),   # no red at 0
#                   (0.5, 1.0, 1.0),   # all channels set to 1.0 at 0.5 to create white
#                   (1.0, 0.8, 0.8)),  # set to 0.8 so its not too bright at 1
#
#         'green': ((0.0, 0.8, 0.8),   # set to 0.8 so its not too bright at 0
#                   (0.5, 1.0, 1.0),   # all channels set to 1.0 at 0.5 to create white
#                   (1.0, 0.0, 0.0)),  # no green at 1
#
#         'blue':  ((0.0, 0.0, 0.0),   # no blue at 0
#                   (0.5, 1.0, 1.0),   # all channels set to 1.0 at 0.5 to create white
#                   (1.0, 0.0, 0.0))   # no blue at 1
#        }
#
# # Create the colormap using the dictionary
# GnRd = colors.LinearSegmentedColormap('GnRd', cdict)
#
# # Make a figure and axes
# fig,ax = plt.subplots(1)
#
# # Some fake data in the range -3 to 3
# dummydata = np.random.rand(4000,18)*6.-3.
#
# # Plot the fake data
# p=ax.pcolormesh(dummydata,cmap=GnRd,vmin=-3,vmax=3)
#
# # Make a colorbar
# fig.colorbar(p,ax=ax)
#
# plt.show()

# import matplotlib.pyplot as plt
# import numpy as np; np.random.seed(0)
# import seaborn as sns; sns.set_theme()
# uniform_data = np.random.rand(4000, 18)
# ax = sns.heatmap(uniform_data, square=True)
# plt.show()

import matplotlib.pyplot as plt
import numpy as np

a = np.random.random((4500, 18))
plt.matshow(a, cmap='hot', interpolation='nearest')
plt.show()
