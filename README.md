# - 四旋翼轨迹跟踪


基于双通道控制机制的轨迹跟踪->运行DisturbanceRejectionControlBasedonDualChannelControlMechanism.m文件
其余控制器的调用都在main.m文件中；通过改变ControllerSelectFlag来选择要运行哪个控制器。attFlag表示姿态旋转，R为旋转矩阵，Q为四元数
各变量名的命名基本符合latex中或word内置unicode的英文表述，需强调变量名后的1d 2d等分别表示对时间的1阶和2阶导数

显示运行结果->运行Display文件

若要在程序运行时查看结果或者要对程序进行调试，请将DisturbanceRejectionControlBasedonDualChannelControlMechanism和main中的全局变量debugFlag置为true

英文版文献
Zhang, Sujie, Xinyu Shi, Xingru Chen and Xiuyu He. “Trajectory Tracking Control Based on Generalized Rodrigues Parameter for Underactuated VTOL UAVs.” Journal of Aerospace Engineering (2024), to be published.

@article{Zhang2024TrajectoryTC,
  title={Trajectory Tracking Control Based on Generalized Rodrigues Parameter for Underactuated VTOL UAVs},
  author={Sujie Zhang and Xinyu Shi and Xingru Chen and Xiuyu He},
  journal={Journal of Aerospace Engineering},
  year={2024},
  url={https://api.semanticscholar.org/CorpusID:272674005}
}
