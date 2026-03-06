from __future__ import annotations
from dataclasses import dataclass
from typing import Optional
import time

from dynamixel_sdk import PortHandler, PacketHandler

def uint32_to_int32(val: int) -> int:
    return val - 0x100000000 if val >= 0x80000000 else val

class DynamxielError(Exception):
    # 统一的Dynamixel错误类型，方便上层捕获
    pass


@dataclass
class DynamixelConfig:
    port:str = 'COM3'  # 替换为你的串口号
    baudrate:int = 57600  # 替换为你的波特率
    dxl_id:int = 1  # 替换为你的Dynamixel ID

    lead_mm_per_rev:float = 2.0  # 线性位移（单位：mm）, 后续精准需要做标定
    gear_ratio:float = 1.67  # 齿轮比
    ticks_per_rev: int = 4096  # 替换为你的Dynamixel的每转动一圈对应的ticks数（单位：ticks）

    max_travel_mm: float = 5.0  # 软限位的最大行程（单位：mm）
    inpos_tol_mm: float = 0.10  # 到位判定的容差（单位：mm）
    timeout_s: float = 10.0  # 等待到位的超时时间（单位：秒）


class XM430Controller:
    #------------------XM430-W210-R控制表地址------------------
    PROTOCOL_VERSION = 2.0
    ADDR_TORQUE_ENABLE = 64
    ADDR_OPERATING_MODE = 11
    ADDR_GOAL_POSITION = 116
    ADDR_PRESENT_POSITION = 132
    MODE_EXTENDED_POSITION = 4

    def __init__(self, config: DynamixelConfig):
        self.cfg = config
        self.port_handler = PortHandler(self.cfg.port)
        self.packet_handler = PacketHandler(self.PROTOCOL_VERSION)
        self.home_ticks: Optional[int] = None
        self._open = False

    def mm_to_ticks(self, delta_mm: float) -> int:
        """将线性位移（mm）转换为Dynamixel的ticks数"""
        revs = delta_mm / self.cfg.lead_mm_per_rev  # 计算转动的圈数
        ticks = revs * self.cfg.ticks_per_rev * self.cfg.gear_ratio  # 转换为ticks数
        return int(round(ticks))
    
    def mm_tol_to_ticks(self, mm_tol: float) -> int:
        """将线性位移的容差（mm）转换为Dynamixel的ticks数"""
        return max(1,abs(self.mm_to_ticks(mm_tol)))
    
    def _write1(self,addr: int, value: int) -> None:
        """向Dynamixel写入一个字节的数据"""
        dxl_comm_result, dxl_error = self.packet_handler.write1ByteTxRx(self.port_handler, self.cfg.dxl_id, addr, value)
        if dxl_comm_result != 0:
            raise DynamxielError(f"comminication error: {self.packet_handler.getTxRxResult(dxl_comm_result)}")
        elif dxl_error != 0:
            raise DynamxielError(f"Dynamixel error: {self.packet_handler.getRxPacketError(dxl_error)}")

    def _write4(self, addr: int, value: int) -> None:
        """向Dynamixel写入四个字节的数据"""
        dxl_comm_result, dxl_error = self.packet_handler.write4ByteTxRx(self.port_handler, self.cfg.dxl_id, addr, value)
        if dxl_comm_result != 0:
            raise DynamxielError(f"comminication error: {self.packet_handler.getTxRxResult(dxl_comm_result)}")
        elif dxl_error != 0:
            raise DynamxielError(f"Dynamixel error: {self.packet_handler.getRxPacketError(dxl_error)}")
    
        
    def _read4(self, addr: int) -> int:
        """从Dynamixel读取四个字节的数据，并返回整数值"""
        dxl_value, dxl_comm_result, dxl_error = self.packet_handler.read4ByteTxRx(self.port_handler, self.cfg.dxl_id, addr)
        if dxl_comm_result != 0:
            raise DynamxielError(f"comminication error: {self.packet_handler.getTxRxResult(dxl_comm_result)}")
        elif dxl_error != 0:
            raise DynamxielError(f"Dynamixel error: {self.packet_handler.getRxPacketError(dxl_error)}")
        return uint32_to_int32(dxl_value)
    
    def open(self)->None:
        """打开串口并初始化Dynamixel"""
        if self._open:
            return
        if not self.port_handler.openPort():
            raise DynamxielError("无法打开串口")
        if not self.port_handler.setBaudRate(self.cfg.baudrate):
            raise DynamxielError("无法设置波特率")
        # 改模式前先Torque Off
        self._write1(self.ADDR_TORQUE_ENABLE, 0)
        self._write1(self.ADDR_OPERATING_MODE, self.MODE_EXTENDED_POSITION)
        self._write1(self.ADDR_TORQUE_ENABLE, 1)
        self._open = True
    
    def close(self)->None:
        """关闭串口"""
        if not self._open:
            return
        try:
            self._write1(self.ADDR_TORQUE_ENABLE, 0)
        except Exception:
            pass
        try:
            self.port_handler.closePort()
        finally:
            self._open = False
    
    def set_home(self)->None:
        # 设置当前位置为HOME，后续move_mm的目标位置都是相对于HOME的增量
        self.home_ticks = self._read4(self.ADDR_PRESENT_POSITION)
        return self.home_ticks
    
    # def clamp_target_to_soft_limits(self, target_ticks: int) -> int:
    #     # 软限位，把target限制在 HOME +/- MAX_TRAVEL_MM 范围内，避免误操作导致机械损伤
    #     if self.home_ticks is None:
    #         return target_ticks  # 没有设置HOME，直接返回原目标
        
    #     max_tricks = abs(self.mm_to_ticks(self.cfg.max_travel_mm))
    #     low = self.home_ticks - max_tricks
    #     high = self.home_ticks + max_tricks
    #     return max(low, min(high, target_ticks))
    
    #------------------- motion------------------------------
    
    def reset_to_home(self, wait=True)->None:
        # 移动到HOME位置
        if self.home_ticks is None:
            raise DynamxielError("HOME position is not set. Cannot reset to home.")
        target = self._validate_target_or_raise(self.home_ticks)
        self._write4(self.ADDR_GOAL_POSITION, target)

        if wait:
            self.wait_until_reached(target)

    def move_from_home_mm(self, mm_from_home: float, wait=True)->None:
        # 从HOME位置移动指定的mm距离
        if self.home_ticks is None:
            raise DynamxielError("HOME position is not set. Cannot move from home.")
        
        target = self.home_ticks + self.mm_to_ticks(mm_from_home)
        target = self._validate_target_or_raise(target)
        self._write4(self.ADDR_GOAL_POSITION, target)
        if wait:
            self._wait_until_reached(target)

    def move_mm_relative(self, delta_mm: float, wait=True)->None:
        # 相对当前位置移动指定的mm距离
        present = self._read4(self.ADDR_PRESENT_POSITION)
        target = present + self.mm_to_ticks(delta_mm)
        target = self._validate_target_or_raise(target) # 应用软限位

        self._write4(self.ADDR_GOAL_POSITION, target)
        if wait:
            self.wait_until_reached(target)
    
    def get_present_ticks(self)->int:
        # 获取当前位置的ticks数
        return self._read4(self.ADDR_PRESENT_POSITION)
    
    def _wait_until_reached(self, target_ticks: int)->None:
        # 等待到位，直到当前位置在目标ticks的threshold范围内，或者超时
        threshold_ticks = self.mm_tol_to_ticks(self.cfg.inpos_tol_mm)

        t0 = time.time()
        while True:
            cur = self._read4(self.ADDR_PRESENT_POSITION)
            if abs(cur - target_ticks) <= threshold_ticks:
                return
            if time.time() - t0 > self.cfg.timeout_s:
                raise DynamxielError(f"Timeout while waiting to reach target. cur={cur}, target={target_ticks}")
            time.sleep(0.02)

    def _validate_target_or_raise(self, target_ticks: int)->int:
        # 验证目标位置是否在软限位范围内，否则抛出异常
        if self.home_ticks is None:
            return target_ticks
        
        max_ticks = abs(self.mm_to_ticks(self.cfg.max_travel_mm))
        low = self.home_ticks - max_ticks
        high = self.home_ticks + max_ticks

        if not (low <= target_ticks <= high):
            raise DynamxielError(f"Target position {target_ticks} ticks is out of soft limits [{low}, {high}] ticks.")    
        return target_ticks