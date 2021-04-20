from energydiagram import ED
import matplotlib.pyplot as plt

diagram = ED()
diagram.add_level(0,'hydroxyacetone')
# 1
diagram.add_level(20.56,'',color='r')
diagram.add_level(23.77,'','l',color='g')
diagram.add_level(12.4,'',color='r')
diagram.add_level(25.16,'','l',color='g')
# 2
diagram.add_level(25.20, '',color='r')
diagram.add_level(31.89, '','l',color='r')
# diagram.add_level(63.81, '','l')
diagram.add_level(28.27, '','l',color='g')
# diagram.add_level(64.36, '','l')
diagram.add_level(3.03, '2-hydroxypropanal',color='r')
diagram.add_level(4.26, '2-hydroxypropanal','l',color='r')
# diagram.add_level(9.07, '3-hydroxypropanal','l')
diagram.add_level(15.11, 'prop-2-ene-1,2-diol','l',color='g')
# diagram.add_level(7.98, '3-hydroxypropanal','l')
#3
diagram.add_level(62.42, '',color='r')
diagram.add_level(59.81, '','l',color='r')
# diagram.add_level(40.62, '','l')
# diagram.add_level(32.85, '','l')
# diagram.add_level(62.72, '','l')
# diagram.add_level(41.05, '','l')
# diagram.add_level(58.76, '','l')
diagram.add_level(65.12, '','l',color='g')
diagram.add_level(19.78, 'acrolein+water',color='r')
diagram.add_level(20.06, 'acrolein+water','l',color='r')
# diagram.add_level(20.76, '','l')
# diagram.add_level(30.90, '','l')
# diagram.add_level(31.70, 'oxetan-2-one+hydrogen','l')
# diagram.add_level(12.31, 'prop-1-ene-1,3-diol','l')
# diagram.add_level(35.95, '','l')
diagram.add_level(39.08, '','l',color='g')
#4
diagram.add_level(71.90, '',color='g')
diagram.add_level(56.49, '','l',color='g')
diagram.add_level(65.12, '','l',color='g')
diagram.add_level(10.26, 'prop-1-ene-1,2-diol',color='g')
diagram.add_level(8.37, 'prop-1-ene-1,2-diol','l',color='g')
diagram.add_level(15.11, 'prop-2-ene-1,2-diol','l',color='g')
#link
diagram.add_link(0,1,color='r')
diagram.add_link(0,2,color='g')
diagram.add_link(1,3,color='r')
diagram.add_link(2,4,color='g')

diagram.add_link(3,6,color='r')
diagram.add_link(3,5,color='r')

diagram.add_link(7,10,color='g')
diagram.add_link(5,8,color='r')
diagram.add_link(6,9,color='r')

diagram.add_link(10,13,color='g')

diagram.add_link(9,12,color='r')
diagram.add_link(8,11,color='r')

diagram.add_link(13,16,color='g')

diagram.add_link(16,17,color='g')
diagram.add_link(12,15,color='r')
diagram.add_link(11,14,color='r')

diagram.add_link(17,20,color='g')
diagram.add_link(4,7,color='g')
diagram.add_link(0,18,color='g')
diagram.add_link(18,21,color='g')
diagram.add_link(16,19,color='g')
diagram.add_link(19,22,color='g')
diagram.plot(show_IDs=False)
plt.show()