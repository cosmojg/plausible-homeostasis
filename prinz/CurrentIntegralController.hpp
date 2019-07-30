// _  _ ____ _    ____ ___ _    
//  \/  |  | |    |  |  |  |    
// _/\_ |__| |___ |__|  |  |___ 
//
// Integral controller, dependent on calcium current
// This controller can control either a synapse
// or a conductance 

#ifndef CURRENTINTEGRALCONTROLLER
#define CURRENTINTEGRALCONTROLLER
#include "mechanism.hpp"
#include <limits>

//inherit controller class spec
class CurrentIntegralController: public mechanism {

protected:
    // flag used to switch between
    // controlling channels and synapses
    // meaning:
    // 0 --- unset, will throw an error
    // 1 --- channels
    // 2 --- synapses
    int control_type = 0;
public:
    // timescales
    double tau_m = std::numeric_limits<double>::infinity();
    double tau_g = 5e3; 

    // mRNA concentration 
    double m = 0;

    // i_Ca
    double tau_filter = std::numeric_limits<double>::quiet_NaN();
    double i_Ca_target = std::numeric_limits<double>::quiet_NaN();
    double i_Ca_error = std::numeric_limits<double>::quiet_NaN();
    double i_Ca_smooth = 0;

    // area of the container this is in
    double container_A;

    // specify parameters + initial conditions for 
    // mechanism that controls a conductance 
    CurrentIntegralController(double tau_m_, double tau_g_, double tau_filter_, double m_, double i_Ca_smooth_, double i_Ca_target_)
    {

        tau_m = tau_m_;
        tau_g = tau_g_;
        tau_filter = tau_filter_;
        m = m_;
        i_Ca_smooth = i_Ca_smooth_;
        i_Ca_target = i_Ca_target_;


        // if (tau_m<=0) {mexErrMsgTxt("[CurrentIntegralController] tau_m must be > 0. Perhaps you meant to set it to Inf?\n");}
        if (tau_g<=0) {mexErrMsgTxt("[CurrentIntegralController] tau_g must be > 0. Perhaps you meant to set it to Inf?\n");}
    }

    
    void integrate(void);

    void checkSolvers(int);

    void connect(conductance *);
    void connect(synapse*);
    void connect(compartment*);

    int getFullStateSize(void);
    int getFullState(double * cont_state, int idx);
    double getState(int);

};


double CurrentIntegralController::getState(int idx)
{
    if (idx == 1) {return m;}
    else if (idx == 2) {return channel->gbar;}
    else {return std::numeric_limits<double>::quiet_NaN();}

}


int CurrentIntegralController::getFullStateSize(){return 2; }


int CurrentIntegralController::getFullState(double *cont_state, int idx)
{
    // give it the current mRNA level
    cont_state[idx] = m;

    idx++;

    // and also output the current gbar of the thing
    // being controller
    if (channel)
    {
      cont_state[idx] = channel->gbar;  
    }
    else if (syn)
    {
        cont_state[idx] = syn->gmax;  
    }
    idx++;
    return idx;
}


void CurrentIntegralController::connect(conductance * channel_)
{

    // connect to a channel
    channel = channel_;


    // make sure the compartment that we are in knows about us
    (channel->container)->addMechanism(this);



    controlling_class = (channel_->getClass()).c_str();

    // attempt to read the area of the container that this
    // controller should be in. 
    container_A  = (channel->container)->A;

    control_type = 1;


}

void CurrentIntegralController::connect(compartment* comp_)
{
    mexErrMsgTxt("[CurrentIntegralController] This mechanism cannot connect to a compartment object");
}

void CurrentIntegralController::connect(synapse* syn_)
{

    // connect to a synpase
    syn = syn_;


    // make sure the compartment that we are in knows about us
    (syn->post_syn)->addMechanism(this);


    // attempt to read the area of the container that this
    // controller should be in. 
    container_A  = (syn->post_syn)->A;

    control_type = 2;

}


void CurrentIntegralController::integrate(void)
{


    switch (control_type)
    {
        case 0:
            mexErrMsgTxt("[CurrentIntegralController] misconfigured controller. Make sure this object is contained by a conductance or synapse object");
            break;


        case 1:

            {
            // if the target is NaN, we will interpret this
            // as the controller being disabled 
            // and do nothing 
            if (isnan(i_Ca_target)) {return;}
          
            // calculate i_Ca_smooth
            i_Ca_smooth += (dt/tau_filter)*((channel->container)->i_Ca_prev-i_Ca_smooth);
            i_Ca_error = i_Ca_smooth - i_Ca_target;

            // integrate mRNA
            m += (dt/tau_m)*(i_Ca_error);

            // mRNA levels below zero don't make any sense
            if (m < 0) {m = 0;}

            // copy the protein levels from this channel
            double gdot = ((dt/tau_g)*(m - channel->gbar*container_A));

            // make sure it doesn't go below zero
            if (channel->gbar_next + gdot < 0) {
                channel->gbar_next = 0;
            } else {
                channel->gbar_next += gdot/container_A;
            }


            break;

            }
        case 2:
            {
            // if the target is NaN, we will interpret this
            // as the controller being disabled 
            // and do nothing 

            if (isnan((syn->post_syn)->Ca_target)) {return;}

            double Ca_error = (syn->post_syn)->Ca_target - (syn->post_syn)->Ca_prev;

            // integrate mRNA
            m += (dt/tau_m)*(Ca_error);

            // mRNA levels below zero don't make any sense
            if (m < 0) {m = 0;}

            // copy the protein levels from this syn
            double gdot = ((dt/tau_g)*(m - syn->gmax*1e-3));

            // make sure it doesn't go below zero
            if (syn->gmax + gdot*1e3 < 0) {
                syn->gmax = 0;
            } else {
                syn->gmax += gdot*1e3;
            }


            break;

            }

        default:
            mexErrMsgTxt("[CurrentIntegralController] misconfigured controller");
            break;

    }


}



void CurrentIntegralController::checkSolvers(int k) {
    if (k == 0){
        return;
    } else {
        mexErrMsgTxt("[CurrentIntegralController] unsupported solver order\n");
    }
}




#endif
