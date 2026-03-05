from dynamixel_sdk import PortHandler, PacketHandler
import time


PORT = 'COM3'  # 替换为你的串口号
BAUDRATE = 57600  # 替换为你的波特率
DXL_ID = 1  # 替换为你的Dynamixel ID

LEAD_MM_PER_REV = 2.0  # 线性位移（单位：mm）, 后续精准需要做标定
GEAR_RATIO = 1.67  # 齿轮比 

#------------------XM430-W210-R控制表地址------------------
PROTOCOL_VERSION = 2.0
ADDR_TORQUE_ENABLE = 64
ADDR_OPERATING_MODE = 11
ADDR_GOAL_POSITION = 116
ADDR_PRESENT_POSITION = 132
MODE_EXTENDED_POSITION = 4

TICKS_PER_REV = 4096  # 替换为你的Dynamixel的每转动一圈对应的ticks数（单位：ticks）

def mm_to_ticks(delta_mm: float) -> int:
    """将线性位移（mm）转换为Dynamixel的ticks数"""
    revs = delta_mm / LEAD_MM_PER_REV  # 计算转动的圈数
    ticks = revs * TICKS_PER_REV * GEAR_RATIO  # 转换为ticks数
    return int(round(ticks))

def write1(packet, port, dxl_id, addr, value):
    """向Dynamixel写入一个字节的数据"""
    dxl_comm_result, dxl_error = packet.write1ByteTxRx(port, dxl_id, addr, value)
    if dxl_comm_result != 0:
        print(f"通信错误: {packet.getTxRxResult(dxl_comm_result)}")
    elif dxl_error != 0:
        print(f"Dynamixel错误: {packet.getRxPacketError(dxl_error)}")


def write4(packet, port, dxl_id, addr, value):
    """向Dynamixel写入四个字节的数据"""
    dxl_comm_result, dxl_error = packet.write4ByteTxRx(port, dxl_id, addr, value)
    if dxl_comm_result != 0:
        print(f"通信错误: {packet.getTxRxResult(dxl_comm_result)}")
    elif dxl_error != 0:
        print(f"Dynamixel错误: {packet.getRxPacketError(dxl_error)}")

def read4(packet, port, dxl_id, addr):
    """从Dynamixel读取四个字节的数据"""
    value, dxl_comm_result, dxl_error = packet.read4ByteTxRx(port, dxl_id, addr)
    if dxl_comm_result != 0:
        print(f"通信错误: {packet.getTxRxResult(dxl_comm_result)}")
        return None
    elif dxl_error != 0:
        print(f"Dynamixel错误: {packet.getRxPacketError(dxl_error)}")
        return None
    return int(value)

def move_mm(packet, port, delta_mm: float, wait=True, threshold_ticks=20, timeout_s=10.0):
    delta_ticks = mm_to_ticks(delta_mm)

    present = read4(packet, port, DXL_ID, ADDR_PRESENT_POSITION)
    target = present + delta_ticks

    write4(packet, port, DXL_ID, ADDR_GOAL_POSITION, target)
    print(f"Move {delta_mm:+.3f} mm -> {delta_ticks:+d} ticks | {present} -> {target}")

    if not wait:
        return

    t0 = time.time()
    while True:
        cur = read4(packet, port, DXL_ID, ADDR_PRESENT_POSITION)
        if abs(cur - target) <= threshold_ticks:
            print(f"Reached: {cur} ticks")
            break
        if time.time() - t0 > timeout_s:
            raise TimeoutError(f"Timeout. cur={cur}, target={target}")
        time.sleep(0.02)


# 主程序
if __name__ == "__main__":
    # 初始化PortHandler和PacketHandler
    port_handler = PortHandler(PORT)
    packet_handler = PacketHandler(PROTOCOL_VERSION)

    # 打开串口
    if not port_handler.openPort():
        raise RuntimeError("无法打开串口")
    # 设置波特率
    if not port_handler.setBaudRate(BAUDRATE):
        raise RuntimeError("无法设置波特率")    

    # 改模式前：Torque Off
    write1(packet_handler, port_handler, DXL_ID, ADDR_TORQUE_ENABLE, 0)
    write1(packet_handler, port_handler, DXL_ID, ADDR_OPERATING_MODE, MODE_EXTENDED_POSITION)
    write1(packet_handler, port_handler, DXL_ID, ADDR_TORQUE_ENABLE, 1)

    try:
        # 示例：移动10mm
        move_mm(packet_handler, port_handler, 0.4)
        time.sleep(2)  # 等待2秒

        # # 示例：移动-5mm
        # move_mm(packet_handler, port_handler, -2)
        # time.sleep(2)  # 等待2秒

    finally:
        # 禁止Dynamixel
        write1(packet_handler, port_handler, DXL_ID, ADDR_TORQUE_ENABLE, 0)
        # 关闭串口
        port_handler.closePort()
        print("Done!")
