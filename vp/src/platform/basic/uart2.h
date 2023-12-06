#ifndef RISCV_ISA_UART2_H
#define RISCV_ISA_UART2_H

#include <tlm_utils/simple_target_socket.h>

#include <systemc>
#include <cstdlib>
#include <queue>
#include <iostream>

#include "core/common/irq_if.h"

template <typename T, int MaxLen, typename Container=std::deque<T>>
    class FixedQueue : public std::queue<T, Container> 
    {
        public:
        void push(const T& value) 
        {
            if (this->size() == MaxLen) {
            this->c.pop_front();
            }
            std::queue<T, Container>::push(value);
        }
    };

struct uart2 : public sc_core::sc_module
{
    tlm_utils::simple_target_socket<uart2> tsock;
    
    interrupt_gateway *plic = 0;
    uint32_t irq_number = 0;
    sc_core::sc_event run_event;

    FixedQueue<char, 8> rx_fifo;
    FixedQueue<char, 8> tx_fifo;
//
    //// memory mapped data frame
    //std::array<char, 8> rx_data;
    //std::array<char, 8> tx_data;

    // memory mapped configuration registers
    u_int32_t tx_data = 0;
    u_int32_t rx_data = 0;
    u_int32_t tx_ctrl = 0;
    u_int32_t rx_ctrl = 0;
    
    std::unordered_map<uint64_t, uint32_t *> addr_to_reg;

    enum {
        TX_DATA_REG_ADDR = 0x00,
        RX_DATA_REG_ADDR = 0x04,
        TX_CTRL_REG_ADDR = 0x08,
        RX_CTRL_REG_ADDR = 0x0C
    };

    SC_HAS_PROCESS(uart2);

    uart2(sc_core::sc_module_name, u_int32_t irq_number) : irq_number(irq_number)
    {
        tsock.register_b_transport(this, &uart2::transport);
        SC_THREAD(run);

        addr_to_reg = {
            {TX_DATA_REG_ADDR, &tx_data},
            {RX_DATA_REG_ADDR, &rx_data},
            {TX_CTRL_REG_ADDR, &tx_ctrl},
            {RX_CTRL_REG_ADDR, &rx_ctrl},
        };
    }

    void transport(tlm::tlm_generic_payload &trans, sc_core::sc_time &delay)
    {
        auto addr = trans.get_address();
		auto cmd = trans.get_command();
		auto len = trans.get_data_length();
		auto ptr = trans.get_data_ptr();

        printf("test");

        (void)delay;
    }



    void run()
    {
        while(true)
        {
            sc_core::wait(sc_core::sc_time(32000, sc_core::SC_NS));
            //create random data
        }
    }
//    {
//        while(true)
//        {
//            run_event.notify(sc_core::sc_time(scaler, sc_core::SC_NS))
//            sc_core::wait(32000, sc_core::SC_NS);
//
//        }
//    }
//
};


#endif //RISCV_ISA_UART2_H