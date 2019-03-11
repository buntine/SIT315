
const int FSR_PIN = 2;
const int LED_PIN = 13;
int LIGHT_PIN = 0;

int photocellReading;
int LEDbrightness;

void pressure_ISR() {
  int pinStatus = digitalRead(FSR_PIN);

  // I'm just tooling for any pressure here in order to enable the LED.
  if (pinStatus > 0) {
    digitalWrite(LED_PIN, HIGH);
    Serial.println("Turn on");
  } else {
    Serial.println("Turn off");
    digitalWrite(LED_PIN, LOW);
  }
}

void setup() 
{
  Serial.begin(9600);
  pinMode(FSR_PIN, INPUT);
  pinMode(13, OUTPUT);

  attachInterrupt(0, pressure_ISR, CHANGE);
}

void loop() 
{
  photocellReading = analogRead(LIGHT_PIN);
 
  Serial.print("Light = ");
  Serial.println(photocellReading);

  // Normalize the photo cell reading over a 0 - 255 range
  // Adapted from: https://learn.adafruit.com/photocells/arduino-code
  photocellReading = 1023 - photocellReading;
  LEDbrightness = map(photocellReading, 0, 1023, 0, 255);

  analogWrite(LED_PIN, LEDbrightness);
  delay(100);
}
