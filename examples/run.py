
import multem

print(multem.is_gpu_available())
print(multem.number_of_gpu_available())

system_conf = multem.SystemConfiguration(device="device", gpu_device=3)
print(system_conf.is_device())
print(system_conf.get_device())

input_multislice = multem.InputMultislice()
input_multislice.set_system_conf(system_conf)
