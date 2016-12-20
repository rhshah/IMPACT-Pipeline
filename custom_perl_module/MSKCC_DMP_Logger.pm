#
# MSKCC_DMP_Log.pm
#
=head1 NAME

MSKCC_DMP_log - loads Log::Log4perl and initializes it 

=head1 SYNOPSIS

 use MSKCC_DMP_Log;

 my $logger = Log::Log4perl->get_logger('My::Class');

and the logger is ready to log, these are the methods in decending
order of verbosity.

 $logger->debug("something pedantic" );
 $logger->info("I did something");
 $logger->warn("I sense danger");
 $logger->error("something is hurting me");
 $logger->fatal("something has killed me");

=head1 DESCRIPTION

MSKCC_DMP_Log is a wrapperish module whose purpose is to initialize 
Log::Log4perl with our standard configuration.  You still get 
loggers from Log4perl.

=head1 VERSION

 Initial Version:     4/30/2013
 Latest Revision:     4/30/2013

=head1 AUTHORS


 Author:              Aijaz Syed
 Contact:             syeda1@mskcc.org

=cut

#####################################################################
# Module Header                                                     #
#####################################################################

package MSKCC_DMP_Logger;

#####################################################################
# Pragma Inclusions                                                 #
#####################################################################

use strict;
use warnings;

#####################################################################
# Declare any other base classes                                    #
#####################################################################
#
# require is used only to generate a less confusing error message 
# when Log::Log4perl isn't available.
#
BEGIN {
    require Log::Log4perl;
};
use Log::Log4perl;
#use base qw(Log::Log4perl);
our @ISA = qw(Log::Log4perl);

#####################################################################
# Module Inclusions                                                 #
#####################################################################
use Carp;
use Cwd;
use File::Path;
use File::Basename;

use MSKCC_DMP_Constants;
use Data::Dumper;
use Log::Log4perl::Layout::PatternLayout;
use Log::Log4perl::Appender;

#####################################################################
# Initialization of module ( main )                                 #
#####################################################################
Log::Log4perl->wrapper_register(__PACKAGE__);

#
# this is called when module loads, at the 'use MSKCC_DMP_Log' statement.
# it allow the optional specification of an alternate configuration
# file.
sub import {
    my ( $self,  $log_type) = @_;

    #if($log_type) {
    $self->_init( '' );
    #}

}

#
# Title      : _init
# Function   : create logging directory if neccessary, init log4perl
# Usage      : _init;
# Args       : none
sub _init {
    my ( $self ) = shift;
    if( ! Log::Log4perl::initialized() ){
        #print "I am here too\n";
        my $configuration;
        if ( get_dev_configuration( \$configuration ) != SUCCESS ){
            carp "can't get configuration to initialize logging!\n";
        }
        Log::Log4perl::init(\$configuration);
        #Log::Log4perl::init("/home/syeda1/Development/Test/log.config");
    }

    #my $logger = Log::Log4perl->get_logger();
    return $self->_adjust_appender('','WARN');
}


sub get_logger {
    my $class = shift;
    my $logger = Log::Log4perl->get_logger(@_);
    my $self = { LOGGER => $logger };
    bless($self,$class);
    return $self;
}


sub fatal { my $self = shift; return $self->{LOGGER}->fatal(@_); }
sub error { my $self = shift; return $self->{LOGGER}->error(@_); }
sub warn  { my $self = shift; return $self->{LOGGER}->warn(@_);  }
sub info  { my $self = shift; return $self->{LOGGER}->info(@_);  }
sub debug { my $self = shift; return $self->{LOGGER}->debug(@_); }
sub trace { my $self = shift; return $self->{LOGGER}->trace(@_); }


#
# Title      : get_prod_configuration
# Function   : creates production configuration for logger
# Usage      : get_prod_configuration <configuration_sReference>;
# Args       : none
sub get_prod_configuration {
    my ($configuration_ref) = @_;
    return get_dev_configuration();
}

#
# Title      : get_dev_configuration
# Function   : creates development configuration for logger
# Usage      : get_dev_configuration <configuration_sReference>;
# Args       : none
sub get_dev_configuration {
    my ($configuration_ref) = @_;
    my $uName = `whoami`;
    chomp($uName);
    my $log_conf = "\
        log4perl.rootLogger              = DEBUG, LOG1\n
        log4perl.appender.LOG1           = Log::Log4perl::Appender::File\n
        log4perl.appender.LOG1.filename  = /home/$uName/IMPACT_central.log\n
        log4perl.appender.LOG1.mode      = append\n
        log4perl.appender.LOG1.layout    = Log::Log4perl::Layout::PatternLayout\n
        log4perl.appender.LOG1.layout.ConversionPattern = %d %p %m %n\n
        log4perl.appender.LOCAL             = Log::Log4perl::Appender::File\n
        log4perl.appender.LOCAL.filename    = Error.log\n
        log4perl.appender.LOCAL.autoflush   = 1\n
        log4perl.appender.LOCAL.mode        = append\n
        log4perl.appender.LOCAL.Threshold   = OFF\n
        log4perl.appender.LOCAL.layout      = Log::Log4perl::Layout::PatternLayout\n
        log4perl.appender.LOCAL.layout.ConversionPattern = %d %p %m %n\n
        log4perl.appender.LOCAL.create_at_logtime = 1\n
    ";
    $$configuration_ref = $log_conf;
    return SUCCESS;
}


sub start_local {
    my ($self, $filePath) = @_;
    my $appender = 'LOCAL';
    
    if (! defined ($filePath) || ! -e $filePath) {
        $filePath = getcwd();
    }

    # Define a file appender
    my $layout = Log::Log4perl::Layout::PatternLayout->new("%d %p %m %n");
    my $file_appender = Log::Log4perl::Appender->new("Log::Log4perl::Appender::File", name => "LOCAL", filename  => $filePath."/IMPACT_local.log");
    $file_appender->layout($layout);
    $self->{LOGGER}->add_appender($file_appender);

    #return $self->get_logger($context);
    my $appenderlist = Log::Log4perl->appenders();
    #print Dumper(%{$appenderlist})."\n";
    if(!defined($appenderlist->{$appender})) {
        #print "Appender list not found...returning failure\n";
        return FAILURE;
    }
    my %TempHash = %Log::Log4perl::Logger::APPENDER_BY_NAME;
    my $app = $TempHash{"LOCAL"};
    my $current_file = $app->filename();
    #my $current_file = "Error.log";
    ##my $filename = $filePath . "/" . $current_file;
    $app->file_switch($current_file);
    my $return_value = $self->_adjust_appender("LOCAL","DEBUG");
    return $return_value;
}


sub _adjust_appender{
    my($self, $appender,$level) = @_;
    my $appenderlist = Log::Log4perl->appenders();
    if(defined($appenderlist->{$appender})) {
        #print "I am here ".$appender."\n";
        my %TempHash = %Log::Log4perl::Logger::APPENDER_BY_NAME;
        my $app = $TempHash{$appender};
        $app->threshold($level);
        my $log = Log::Log4perl->get_logger('');
        my $root_level = $log->level();
    } else {
        if( ! Log::Log4perl::initialized() ){
            carp "Error initializing the local logger..  ";
            return FAILURE;
        }
    }
    return SUCCESS;
}


sub stop_local {
    my $self = shift;
    my $return_value;
    # Turn on the optional appender
    $return_value = $self->_adjust_appender('Optional','OFF');
    return $return_value;
}

 
#
# A seemingly gratuitous statement, included so that the interpreter
# doesn't complain. return true (1) value so perl knows module loaded succesfully
1;
