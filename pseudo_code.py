# Structure of code for electrmechanical stretcher
#
# Modules:
#     Python:
#         Arduino serial reader (load cell + stepper state)
#         Ohmmeter serial reader
#         Overall data collector and file writer
#     C++:
#         Loadcell reader library
#         Stepper driver
#         Linear actuator
#         Other sensors?


Process flow:
    C++(Arduino):
        1) Initialise
            a) Stepper - step size, ustepping
            b) Loadcell -
            c) Serial layer - Optimise speed?, send/recieve protocol setup
        2) Process
            a) Read Serial
                -> Run stepper profile
                -> Read loadcell val
                -> Read stepper pos
                (-> timestamp each message?)
            b) Write serial
                -> Send loadcell values
                -> Send stepper pos
                -> Send time
                -> Send error msgs
            c) Check stepper profile state (timer interupt, on a 'per step' basis)

    Python:
        1) Initialise
            a) Linear actuator - Max force, total travel, thread pitch, stepper vals, cmd interface
            b) Loadcell - conversion factors
            c) Serial interfaces - Arduino + Ohmmeter
        2) Process
            a) Setup experiment sequence - velocity profile, setup ohmmeter params.
            b) Loop until sequence finished
                i) Read arduino time (to compare with PC time, get time when serial packet sent and recieved to estimate send time offset)
                ii) Read loadcell + time
                iii) Read stepper pos + time
                iv) Read Ohmmeter val + time
                v) Read loadcell + time? not req?
                vi) Read stepper pos + time? not req?
                vii) if any errors occur pause/stop/react
            c) Write data into a .csv file in time order
