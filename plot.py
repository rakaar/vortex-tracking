vortices = [
     [          65 ,          32 ],
 [          65 ,         160 ],
 [         193 ,          32 ],
 [         193 ,         160 ]
]
x_cordinates = [x[0] for x in vortices]
y_cordinates = [x[1] for x in vortices]

import matplotlib.pyplot as plt
plt.scatter(x_cordinates, y_cordinates, )
plt.show()
